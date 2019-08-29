# Mouse_fmriPype: Automated BOLD preprocessing pipeline for rodent fMRI analysis.

![Processing Schema](https://github.com/Gab-D-G/pics/blob/master/processing_schema.jpg)

## Input data folder structure
* input_folder/subject_id/ses-#/bold/subject_id_ses-#_run-#_bold.nii.gz
* input_folder/subject_id/ses-#/anat/subject_id_ses-#_anat.nii.gz
* input_folder/data_info.csv (with columns header: group,subject_id,num_session,num_run)


<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
 <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
 <meta name="Author" content="Made by 'tree'">
 <meta name="GENERATOR" content="$Version: $ tree v1.6.0 (c) 1996 - 2011 by Steve Baker, Thomas Moore, Francesc Rocher, Kyosuke Tokoro $">
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
	<h1>Example Directory Tree</h1><p>
	<a href="nii_inputs/">nii_inputs/</a><br>
	└── <a href="nii_inputs//nii_inputs/">nii_inputs</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/data_info.csv">data_info.csv</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/">jgrAesAWc11L</a><br>
	&nbsp;&nbsp;&nbsp; │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; │   │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/anat/jgrAesAWc11L_ses-1_anat.nii.gz">jgrAesAWc11L_ses-1_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/bold/jgrAesAWc11L_ses-1_run-1_bold.nii.gz">jgrAesAWc11L_ses-1_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/bold/jgrAesAWc11L_ses-1_run-2_bold.nii.gz">jgrAesAWc11L_ses-1_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-1/bold/jgrAesAWc11L_ses-1_run-3_bold.nii.gz">jgrAesAWc11L_ses-1_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/">ses-2</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/anat/jgrAesAWc11L_ses-2_anat.nii.gz">jgrAesAWc11L_ses-2_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/bold/jgrAesAWc11L_ses-2_run-1_bold.nii.gz">jgrAesAWc11L_ses-2_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/bold/jgrAesAWc11L_ses-2_run-2_bold.nii.gz">jgrAesAWc11L_ses-2_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11L/ses-2/bold/jgrAesAWc11L_ses-2_run-3_bold.nii.gz">jgrAesAWc11L_ses-2_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/">jgrAesAWc11R</a><br>
	&nbsp;&nbsp;&nbsp; │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; │   │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/anat/jgrAesAWc11R_ses-1_anat.nii.gz">jgrAesAWc11R_ses-1_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/bold/jgrAesAWc11R_ses-1_run-1_bold.nii.gz">jgrAesAWc11R_ses-1_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/bold/jgrAesAWc11R_ses-1_run-2_bold.nii.gz">jgrAesAWc11R_ses-1_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-1/bold/jgrAesAWc11R_ses-1_run-3_bold.nii.gz">jgrAesAWc11R_ses-1_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/">ses-2</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/anat/jgrAesAWc11R_ses-2_anat.nii.gz">jgrAesAWc11R_ses-2_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/bold/jgrAesAWc11R_ses-2_run-1_bold.nii.gz">jgrAesAWc11R_ses-2_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/bold/jgrAesAWc11R_ses-2_run-2_bold.nii.gz">jgrAesAWc11R_ses-2_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R/ses-2/bold/jgrAesAWc11R_ses-2_run-3_bold.nii.gz">jgrAesAWc11R_ses-2_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/">jgrAesAWc11R1L</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/">ses-1</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/anat/jgrAesAWc11R1L_ses-1_anat.nii.gz">jgrAesAWc11R1L_ses-1_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/bold/jgrAesAWc11R1L_ses-1_run-1_bold.nii.gz">jgrAesAWc11R1L_ses-1_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/bold/jgrAesAWc11R1L_ses-1_run-2_bold.nii.gz">jgrAesAWc11R1L_ses-1_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-1/bold/jgrAesAWc11R1L_ses-1_run-3_bold.nii.gz">jgrAesAWc11R1L_ses-1_run-3_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/">ses-2</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/anat/">anat</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; │   └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/anat/jgrAesAWc11R1L_ses-2_anat.nii.gz">jgrAesAWc11R1L_ses-2_anat.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/bold/">bold</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/bold/jgrAesAWc11R1L_ses-2_run-1_bold.nii.gz">jgrAesAWc11R1L_ses-2_run-1_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; ├── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/bold/jgrAesAWc11R1L_ses-2_run-2_bold.nii.gz">jgrAesAWc11R1L_ses-2_run-2_bold.nii.gz</a><br>
	&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; └── <a href="nii_inputs//nii_inputs/jgrAesAWc11R1L/ses-2/bold/jgrAesAWc11R1L_ses-2_run-3_bold.nii.gz">jgrAesAWc11R1L_ses-2_run-3_bold.nii.gz</a><br>
	<br><br>
	</p>
	<p>

22 directories, 25 files
	<br><br>
	</p>
	<hr>
	<p class="VERSION">
		 tree v1.6.0 © 1996 - 2011 by Steve Baker and Thomas Moore <br>
		 HTML output hacked and copyleft © 1998 by Francesc Rocher <br>
		 Charsets / OS/2 support © 2001 by Kyosuke Tokoro
	</p>
</body>
</html>


\nFor questions or interests in using the pipeline, contact gabriel.desrosiers-gregoire@mail.mcgill.ca

## References

* fmriprep - https://github.com/poldracklab/fmriprep
* minc-toolkit v2 - https://github.com/BIC-MNI/minc-toolkit-v2
* minc-stuffs - https://github.com/Mouse-Imaging-Centre/minc-stuffs
* minc2-simple - https://github.com/vfonov/minc2-simple
* pydpiper - https://github.com/Mouse-Imaging-Centre/pydpiper
