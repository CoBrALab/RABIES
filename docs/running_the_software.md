# Running The Software

## Input data in BIDS standard

The input dataset must be organized according to the [BIDS data structure](https://bids.neuroimaging.io/){cite}`Gorgolewski2016-zm`. RABIES will iterate through all subjects found to contain a functional file (see section on BIDS filters below), and will also iterate according to sessions and runs found within each subject if available. If anatomical scans are used for preprocessing (i.e. not using `--bold_only`), each functional scan will be matched to one corresponding anatomical scan of the same subject/session (using BIDS filters for the anatomical image, see below).

### Directory structure for an example dataset
* Our [example dataset](http://doi.org/10.5281/zenodo.8349029) has the following BIDS structure: 

<!DOCTYPE html>
<html>
<head>
 <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
 <meta name="Author" content="Made by 'tree'">
 <meta name="GENERATOR" content="$Version: $ tree v1.7.0 (c) 1996 - 2014 by Steve Baker, Thomas Moore, Francesc Rocher, Florian Sesser, Kyosuke Tokoro $">
  <!-- 
  BODY { font-family : ariel, monospace, sans-serif; }
  P { font-weight: normal; font-family : ariel, monospace, sans-serif; color: black; background-color: transparent;}
  B { font-weight: normal; color: black; background-color: transparent;}
  A:visited { font-weight : normal; text-decoration : none; background-color : transparent; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }
  A:link    { font-weight : normal; text-decoration : none; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }
  A:hover   { color : #000000; font-weight : normal; text-decoration : underline; background-color : yellow; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }
  A:active  { color : #000000; font-weight: normal; background-color : transparent; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }
  .VERSION { font-size: small; font-family : arial, sans-serif; }
  .NORM  { color: black;  background-color: transparent;}
  .FIFO  { color: purple; background-color: transparent;}
  .CHAR  { color: yellow; background-color: transparent;}
  .DIR   { color: blue;   background-color: transparent;}
  .BLOCK { color: yellow; background-color: transparent;}
  .LINK  { color: aqua;   background-color: transparent;}
  .SOCK  { color: fuchsia;background-color: transparent;}
  .EXEC  { color: green;  background-color: transparent;}
  -->
</head>
<body>
        <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip">https://zenodo.org/record/8349029/preview/test_dataset.zip</a><br>
        ├── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG001/">sub-PHG001</a><br>
        │   └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG001/ses-3/">ses-3</a><br>
        │   &nbsp;&nbsp;&nbsp; ├── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG001/ses-3/anat/">anat</a><br>
        │   &nbsp;&nbsp;&nbsp; │   ├── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG001/ses-3/anat/sub-PHG001_ses-3_acq-RARE_T2w.json">sub-PHG001_ses-3_acq-RARE_T2w.json</a><br>
        │   &nbsp;&nbsp;&nbsp; │   └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG001/ses-3/anat/sub-PHG001_ses-3_acq-RARE_T2w.nii.gz">sub-PHG001_ses-3_acq-RARE_T2w.nii.gz</a><br>
        │   &nbsp;&nbsp;&nbsp; └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG001/ses-3/func/">func</a><br>
        │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG001/ses-3/func/sub-PHG001_ses-3_task-rest_acq-EPI_run-1_bold.json">sub-PHG001_ses-3_task-rest_acq-EPI_run-1_bold.json</a><br>
        │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG001/ses-3/func/sub-PHG001_ses-3_task-rest_acq-EPI_run-1_bold.nii.gz">sub-PHG001_ses-3_task-rest_acq-EPI_run-1_bold.nii.gz</a><br>
        └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG002/">sub-PHG002</a><br>
        &nbsp;&nbsp;&nbsp; └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG002/ses-3/">ses-3</a><br>
        &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG002/ses-3/anat/">anat</a><br>
        &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   ├── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG002/ses-3/anat/sub-PHG002_ses-3_acq-RARE_T2w.json">sub-PHG002_ses-3_acq-RARE_T2w.json</a><br>
        &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG002/ses-3/anat/sub-PHG002_ses-3_acq-RARE_T2w.nii.gz">sub-PHG002_ses-3_acq-RARE_T2w.nii.gz</a><br>
        &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG002/ses-3/func/">func</a><br>
        &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG002/ses-3/func/sub-PHG002_ses-3_task-rest_acq-EPI_run-1_bold.json">sub-PHG002_ses-3_task-rest_acq-EPI_run-1_bold.json</a><br>
        &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="https://zenodo.org/record/8349029/preview/test_dataset.zip/sub-PHG002/ses-3/func/sub-PHG002_ses-3_task-rest_acq-EPI_run-1_bold.nii.gz">sub-PHG002_ses-3_task-rest_acq-EPI_run-1_bold.nii.gz</a><br>
        <br><br>
        </p>
        <p>

8 directories, 8 files
        <br><br>
        </p>
        <hr>
</body>
</html>

### BIDS filters to identify functional and structural images
By default, RABIES will use the 'bold' or 'cbv' suffix to identify functional scans and the 'T1w' or 'T2w' suffix for structural scans. Files which don't match the BIDS filters are ignored. However, the BIDS filters can also be customized with the `--bids_filter` parameter during the preprocessing stage. This can be useful for instance if the default is not enough to find the right set of scans. The custom BIDS filter must be formated into a JSON file with the functional filter under 'func' and structural filter under 'anat' (see example below for the default parameters):
```json
{
    "func": {
        "suffix":["bold","cbv"]
    },
    "anat": {
        "suffix":["T1w","T2w"]
    }
}
```

## Command Line Interface

RABIES is executed using a command line interface, within a terminal. The software is divided into three main processing stages: preprocessing, confound correction and analysis. Accordingly, the command line interface allows for three different mode of execution, corresponding to the processing stages. So first, when executing the software, one of the processing stage must be selected. Below you can find the general --help message printed with `rabies --help`, which provides a summary of each processing stage together with options for parallel processing and memory management. Then, the --help associated to each processing stage, i.e. `preprocess`, `confound_correction` and `analysis`, describes in more detail the various parameters available to adapt image processing according to the user needs. Click on the corresponding --help to expand:

<details><summary><b>rabies --help</b></summary>
<p>

```{program-output} rabies --help
```

</p>
</details>

<details><summary><b>rabies preprocess --help</b></summary>
<p>

```{program-output} rabies preprocess --help
```

</p>
</details>

<details><summary><b>rabies confound_correction --help</b></summary>
<p>

```{program-output} rabies confound_correction --help
```

</p>
</details>


<details><summary><b>rabies analysis --help</b></summary>
<p>

```{program-output} rabies analysis --help
```

</p>
</details>


## Example execution syntax
The following section provides examples describing the basic syntax for running the RABIES command line interface.


**preprocess**
```sh
rabies -p MultiProc preprocess input_BIDS/ preprocess_outputs/ --apply_STC --TR 1.2 --commonspace_reg masking=true,brain_extraction=false,template_registration=SyN,fast_commonspace=false
```
First, we have to preprocess the dataset before it can be analyzed. In this example, we are running the RABIES preprocessing on the dataset found in the `input_BIDS/` folder, formatted according to the BIDS standard, and the outputs from RABIES are stored in the `preprocess_outputs/` folder. Additional execution parameters were specified: 
* `-p MultiProc` will execute the pipeline in parallel using the local threads available. Notice that this parameter is specified before the processing stage, because it is one of the `Execution Options` affiliated to the `rabies --help`.
* `--apply_STC` is a boolean variable which, when selected, will apply slice timing correction during preprocessing, which is not applied by default in RABIES.
* `--TR 1.2` specifies the repetition time (TR) of the fMRI images that are processed, which must be defined to apply slice timing correction appropriately. Notice that this parameter must be provided with an argument, here `1.2` for TR = 1.2sec, and this is done by writing down the argument with a space dividing the associated parameter.
* `--commonspace_reg masking=true,brain_extraction=false,template_registration=SyN,fast_commonspace=false` this argument manages the options for the commonspace registration step. Some arguments, including `--commonspace_reg`, take multiple parameters as input, where each parameter-value pairs follow the syntax of `parameter=value`. In this case, we are using the masking optino with `masking=true`, which use available brain masks from inhomogeneity correction to drive the registration operations, and we specify a non-linear registration to the commonspace template with `template_registration=SyN`.

**confound_correction**
```sh
rabies -p MultiProc confound_correction preprocess_outputs/ confound_correction_outputs/ --conf_list WM_signal CSF_signal vascular_signal mot_6 --smoothing_filter 0.3 
```
Next, after completing preprocessing, in most cases the data should be corrected for potential confounds prior to analysis. This is done in the confound correction stage, where confounds are modelled and regressed from the data. In this example we correct the preprocessed data found in the `preprocess_outputs/` folder and store the cleaned outputs in the `confound_correction_outputs/` folder. Among the range of options available for confound correction, we define in this example three parameters:
* `--conf_list` is the option to regress nuisance timeseries from the data, i.e., confound regression. This parameter takes a list as input, where each argument in the list is seperated by a space as follow `WM_signal CSF_signal mot_6`. This list defines which nuisance timeseries are going to model confounds during confound regression, in this case, the WM and CSF mean signals together with the 6 rigid realignment parameters from head motion realignment.
* `--smoothing_filter` will additionally apply Gaussian spatial smoothing, where in this case, a filter size of `0.3` mm is specified.

**analysis**
```sh
rabies -p MultiProc analysis confound_correction_outputs analysis_outputs/ --group_ica apply=true,dim=30,random_seed=1
```
Finally, after conducting preprocessing and confound correction, certain analyses can be run within RABIES. In this case, the cleaned outputs found in `confound_correction_outputs/` are going to be analyzed, with analysis outputs found in `analysis_outputs/`. We perform a group independent component analysis (ICA) with 30 components by providing `--group_ica apply=true,dim=30,random_seed=1` to the command.

## Execution syntax with containerized installation (Singularity and Docker)

Containers are independent computing environments which have their own dependencies installed to ensure consistent and reliable
execution of the software across computing platforms. Singularity containers, as opposed to Docker, can be exported to remote high-performance computing platforms (e.g. computecanada). The main difference in execution syntax when running a container, as opposed to the examples above, is that the paths between the local environment where the data is stored must be 'linked' to the container's internal paths. All relevant directories containing data that will be used by RABIES must be related to a container internal path, and this is done using `-B` for Singularity and `-v` for Docker. See below for examples:

### Singularity execution

**preprocess**
```sh
singularity run -B $PWD/input_BIDS:/input_BIDS:ro \
-B $PWD/preprocess_outputs:/preprocess_outputs/ \
/path_to_singularity_image/rabies.sif -p MultiProc preprocess /input_BIDS/ /preprocess_outputs/ --apply_STC --TR 1.2 --commonspace_reg masking=true,brain_extraction=false,template_registration=SyN,fast_commonspace=false
```
Singularity containers are stored in image files, for instance `rabies.sif`. `singularity run /path_to_singularity_image/rabies.sif` will execute the image, in this case the RABIES pipeline, and the same rules for the command line interface then apply as previously demonstrated. However, the container must gain access to the relevant folders for running RABIES, in this case an input folder and an output folder, and this is done with `-B`:
* `-B $PWD/input_BIDS:/input_BIDS:ro`: this argument relates the BIDS input folder found in `$PWD/input_BIDS` to an internal path to the container, which we call `/input_BIDS`. The inputs are thus accessed according to this path in the RABIES arguments with `/input_BIDS/`. the `:ro` means that the container is only provided reading permissions at this location.
* `-B $PWD/preprocess_outputs:/preprocess_outputs/`: same as with the `/input_BIDS/`, but now we are relating a desired output directory `$PWD/preprocess_outputs` to `/preprocess_outputs`, and the container has writing permissions at this path since `:ro` is not present.


**confound_correction**
```sh
singularity run -B $PWD/input_BIDS:/input_BIDS:ro \
-B $PWD/preprocess_outputs:/preprocess_outputs/ \
-B $PWD/confound_correction_outputs:/confound_correction_outputs/ \
/path_to_singularity_image/rabies.sif -p MultiProc confound_correction /preprocess_outputs/ /confound_correction_outputs/ --conf_list WM_signal CSF_signal vascular_signal mot_6 --smoothing_filter 0.3 
```
The required paths are similarly provided for the confound correction stage. Note here that the path to `$PWD/input_BIDS` is still linked to the container, even though it is not explicitely part of the arguments during the confound correction call. This is necessary since the paths used in the preprocessing steps still need to be accessed at later stages, and there will be an error if the paths are not kept consistent across processing steps.

**analysis**
```sh
singularity run -B $PWD/input_BIDS:/input_BIDS:ro \
-B $PWD/preprocess_outputs:/preprocess_outputs/ \
-B $PWD/confound_correction_outputs:/confound_correction_outputs/ \
-B $PWD/analysis_outputs:/analysis_outputs/ \
/path_to_singularity_image/rabies.sif -p MultiProc analysis /confound_correction_outputs /analysis_outputs/ --group_ica apply=true,dim=30,random_seed=1
```
The same logic applies at the analysis stage.
<br/>

### Docker execution
```sh
docker run -it --rm --user $(id -u) \
-v $PWD/input_BIDS:/input_BIDS:ro \
-v $PWD/preprocess_outputs:/preprocess_outputs/ \
gabdesgreg/rabies:tagname -p MultiProc preprocess /input_BIDS/ /preprocess_outputs/ --apply_STC --TR 1.2 --commonspace_reg masking=true,brain_extraction=false,template_registration=SyN,fast_commonspace=false
```
The syntax in Docker is very similar to Singularity, except that `-B` is replaced by `-v`, and further parameters may be needed (e.g. `-it`, `--rm`). `--user $(id -u)` can be added to mitigate writing permission issues when using Docker. Note that 'tagname' should be replaced by the proper RABIES version you are using (e.g. 0.4.8).


## Additional Resources

* **Workshop and tutorial for RABIES:**
    * [Hands-on tutorial](https://github.com/CoBrALab/RABIES_tutorial) on RABIES (originally developed for and presented at the INCF Neuroinformatics Assembly 2023).
    * A workshop providing a complete software overview was [recorded and posted online](https://www.youtube.com/watch?v=LZohKlUgycc&t=2766s&ab_channel=DouglasResearchCentre) on February 2023.
* Conversion from Bruker raw to Nifti formats can be handled with [BrkRaw](https://brkraw.github.io/) (consult [associated documentation](https://github.com/CoBrALab/documentation/wiki/bruker2nifti-conversion) from the CoBrALab)
* [CoBrALab recommendations](https://github.com/CoBrALab/documentation/wiki/Running-RABIES-on-niagara) for using compute canada.

