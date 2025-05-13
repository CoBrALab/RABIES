./docs/conf.py - some path stuff but nothing that requires changing
./minc-toolkit-extras/ants_generate_iterations.py - nothing
./minc-toolkit-extras/t2star_fit.py - nothing
./minc-toolkit-extras/t2star_fit_simpleitk.py - nothing
./optimized_antsMultivariateTemplateConstruction/modelbuild_averager.py - has os.path but is a submodule
./rabies/analysis_pkg/analysis_functions.py - already had some pathlib, had lots of os.path.abspath, changed to pathlib.Path(path).absolute()
./rabies/analysis_pkg/analysis_math.py - nothing
./rabies/analysis_pkg/analysis_wf.py - few, minor, using is_file() method now too
./rabies/analysis_pkg/diagnosis_pkg/analysis_QC.py - nothing
./rabies/analysis_pkg/diagnosis_pkg/diagnosis_functions.py - done
./rabies/analysis_pkg/diagnosis_pkg/diagnosis_wf.py - basename -> .name
./rabies/analysis_pkg/diagnosis_pkg/**init**.py - nothing
./rabies/analysis_pkg/diagnosis_pkg/interfaces.py - figure_path needed to be converted back to a string to be used by `savefig` function, this should be checked
also using pathlib `mkdir` instead of `os.mkdir`
around line 335, analysis_QC_network there's some path stuff that may need to be fixed
./rabies/analysis_pkg/**init**.py - done
./rabies/analysis_pkg/main_wf.py - done
./rabies/analysis_pkg/utils.py - done

###

./rabies/boilerplate.py - nothing
./rabies/confound_correction_pkg/confound_correction.py - replaced os.path with pathlib. there are a lot of strings as file paths. may be an issue
./rabies/confound_correction_pkg/**init**.py - nothing
./rabies/confound_correction_pkg/main_wf.py - is_dir() used, figure path to a string again. re-visit, may need to be re-done

###

./rabies/confound_correction_pkg/mod_ICA_AROMA/classification_plots.py replaced os.path.join with pathlib `/` syntax
./rabies/confound_correction_pkg/mod_ICA_AROMA/ICA_AROMA_functions.py changed to use `with open`, remvoed txt.close(), files will close correctly if error occurs
currently uses camel case insteead of snake case, this shold be changed
a lot of changes in this file, review required

./rabies/confound_correction_pkg/mod_ICA_AROMA/ICA_AROMA.py
./rabies/confound_correction_pkg/mod_ICA_AROMA/ica-aroma-via-docker.py
./rabies/confound_correction_pkg/mod_ICA_AROMA/**init**.py
./rabies/confound_correction_pkg/utils.py
./rabies/**init**.py
./rabies/parser.py
./rabies/preprocess_pkg/bold_main_wf.py
./rabies/preprocess_pkg/bold_ref.py
./rabies/preprocess_pkg/commonspace_reg.py
./rabies/preprocess_pkg/hmc.py
./rabies/preprocess_pkg/inho_correction.py
./rabies/preprocess_pkg/**init**.py
./rabies/preprocess_pkg/main_wf.py
./rabies/preprocess_pkg/preprocess_visual_QC.py
./rabies/preprocess_pkg/registration.py
./rabies/preprocess_pkg/resampling.py
./rabies/preprocess_pkg/stc.py
./rabies/preprocess_pkg/utils.py
./rabies/run_main.py
./rabies/utils.py
./rabies/**version**.py
./rabies/visualization.py
./scripts/debug_workflow.py
./scripts/error_check_rabies.py
./scripts/gen_DSURQE_masks.py
./scripts/preprocess_scripts/multistage_otsu_cor.py
./setup.py
