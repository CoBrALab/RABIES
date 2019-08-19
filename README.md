Automated BOLD preprocessing pipeline for rodent fMRI analysis.

![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/processing_schema.jpg)

Input data structure: input_folder/subject_id/ses-#/bold/subject_id_ses-#_run-#_bold.nii.gz, input_folder/subject_id/ses-#/anat/subject_id_ses-#_anat.nii.gz,
                      input_folder/data_info.csv
data_info.csv: columns group,subject_id,num_session,num_run

![Input Data Structure](https://github.com/Gab-D-G/pics/blob/master/input_structure.jpg)

For questions or interests in using the pipeline, contact gabriel.desrosiers-gregoire@mail.mcgill.ca


**References**
fmriprep - https://github.com/poldracklab/fmriprep
minc-toolkit v2 - https://github.com/BIC-MNI/minc-toolkit-v2
minc-stuffs - https://github.com/Mouse-Imaging-Centre/minc-stuffs
minc2-simple - https://github.com/vfonov/minc2-simple
pydpiper - https://github.com/Mouse-Imaging-Centre/pydpiper
