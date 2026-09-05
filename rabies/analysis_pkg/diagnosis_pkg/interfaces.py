import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
from rabies.analysis_pkg.diagnosis_pkg import diagnosis_functions
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


class ScanDiagnosisInputSpec(BaseInterfaceInputSpec):
    cleaned_bold_file = File(exists=True, mandatory=True,
        desc="Cleaned/confound-corrected 4D EPI timeseries for this scan, in the analysis space (native or commonspace).")
    name_source = File(exists=True, mandatory=True,
        desc="The original raw EPI file, used only to name outputs consistently with the rest of RABIES.")
    CR_data_dict = traits.Dict(
        desc="Confound correction QC dict for this scan (FD_trace, tDOF, frame_mask, TR, ...); 'voxelwise_mean' is also read from it if remove_EPI_avg is True.")
    remove_EPI_avg = traits.Bool(False, usedefault=True,
        desc="Whether the average was kept in the cleaned timeseries and must be removed prior to analysis.")
    VE_file = File(exists=True, mandatory=True, desc="Variance explained (R^2) map from confound regression.")
    STD_file = File(exists=True, mandatory=True, desc="Temporal standard deviation map on the cleaned timeseries.")
    CR_STD_file = File(exists=True, mandatory=True, desc="Temporal standard deviation map on the predicted confound timeseries.")

    # subject-space definition (native or commonspace, wherever confound correction ran)
    sub_mask_file = File(exists=True, mandatory=True, desc="Brain mask in the subject's analysis space.")
    sub_WM_mask_file = traits.Either(File(exists=True), None, desc="WM mask in the subject's analysis space.")
    sub_CSF_mask_file = traits.Either(File(exists=True), None, desc="CSF mask in the subject's analysis space.")

    # commonspace definition (always commonspace, used for display and group comparison)
    common_mask_file = File(exists=True, mandatory=True, desc="Brain mask in commonspace.")
    common_anat_ref_file = File(exists=True, mandatory=True, desc="Anatomical template in commonspace.")

    analysis_files_dict = traits.Dict(
        desc="A dictionary regrouping relevant output files from analysis.")
    prior_bold_idx = traits.List(
        desc="The index for the ICA components that correspond to bold sources.")
    prior_confound_idx = traits.List(
        desc="The index for the ICA components that correspond to confounding sources.")
    plot_seed_frequencies = traits.Dict(
        desc="Dictionary that pairs a seed name and its index within --seed_list to plot the frequency spectrum")
    figure_format = traits.Str(
        desc="Select file format for figures.")
    nativespace_analysis = traits.Bool(
        desc="Whether input timeseries are in nativespace.")
    native_to_common_transforms = traits.List(
        desc="Transforms to move nativespace computations to commonspace")
    native_to_common_inverses = traits.List(
        desc="Inverses to move nativespace computations to commonspace")
    interpolation = traits.Int(
        desc="Provide an SITK interpolator.")
    brainmap_percent_threshold = traits.Float(
        desc="Input percentage value for thresholding images.")
    rabies_data_type = traits.Int(
        desc="Integer specifying SimpleITK data type.")


class ScanDiagnosisOutputSpec(TraitedSpec):
    figure_temporal_diagnosis = File(
        exists=True, desc="Output figure from the scan diagnosis")
    figure_spatial_diagnosis = File(
        exists=True, desc="Output figure from the scan diagnosis")
    temporal_info = traits.Dict(
        desc="A dictionary regrouping the temporal features.")
    spatial_info = traits.Dict(
        desc="A dictionary regrouping the spatial features.")


class ScanDiagnosis(BaseInterface):
    """
    Extracts several spatial and temporal features on the target scan.
    Spatial features include tSTD, CRsd (variance fitted from confound regression), CR-R^2,
    and network maps from specified BOLD priors at the indices of prior_bold_idx.
    Temporal features include grayplot, 6 motion parameters, framewise displacement,
    DVARS, WM/CSV/edge mask timecourses, CR-R^2, and the average amplitude of BOLD and
    confound components seperately.
    """

    input_spec = ScanDiagnosisInputSpec
    output_spec = ScanDiagnosisOutputSpec

    def _run_interface(self, runtime):
        import pathlib
        from rabies.analysis_pkg.utils import load_resample_analysis_maps

        figure_format = self.inputs.figure_format
        CR_data_dict = self.inputs.CR_data_dict
        nativespace_analysis = self.inputs.nativespace_analysis
        filename_split = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")[0]

        # subject-space data (matches wherever confound correction ran): timeseries + WM/CSF/edge indices.
        # no resampling is needed here -- native_WM_mask/native_brain_mask etc. are already in this space --
        # so anat_ref_file is passed as sub_mask_file itself, purely as a placeholder.
        sub_space_maps_loaded = load_resample_analysis_maps(
            self.inputs.sub_mask_file, self.inputs.sub_mask_file,
            cleaned_bold_file=self.inputs.cleaned_bold_file, CR_data_dict=CR_data_dict, remove_EPI_avg=self.inputs.remove_EPI_avg,
            WM_mask_file=self.inputs.sub_WM_mask_file, CSF_mask_file=self.inputs.sub_CSF_mask_file, compute_edge_idx=True,
            output_prefix=filename_split)

        VE_spatial = sitk.GetArrayFromImage(sitk.ReadImage(self.inputs.VE_file))[sub_space_maps_loaded['volume_indices']]
        temporal_std = sitk.GetArrayFromImage(sitk.ReadImage(self.inputs.STD_file))[sub_space_maps_loaded['volume_indices']]
        predicted_std = sitk.GetArrayFromImage(sitk.ReadImage(self.inputs.CR_STD_file))[sub_space_maps_loaded['volume_indices']]

        # convert to an integer list
        prior_bold_idx = [int(i) for i in self.inputs.prior_bold_idx]
        prior_confound_idx = [int(i) for i in self.inputs.prior_confound_idx]

        if nativespace_analysis:
            resampling_specs = {'transforms':self.inputs.native_to_common_transforms,
                                'inverses':self.inputs.native_to_common_inverses,
                                'interpolation':self.inputs.interpolation,
                                'rabies_data_type':self.inputs.rabies_data_type,
                                }
        else:
            resampling_specs = {}

        temporal_info, spatial_info = diagnosis_functions.compute_spatiotemporal_features(
            self.inputs.name_source, CR_data_dict, VE_spatial, temporal_std, predicted_std,
            self.inputs.sub_mask_file, sub_space_maps_loaded, self.inputs.common_mask_file, self.inputs.common_anat_ref_file,
            self.inputs.analysis_files_dict, prior_bold_idx, prior_confound_idx,
            nativespace_analysis=nativespace_analysis, resampling_specs=resampling_specs)

        temporal_fig_list, fig2 = diagnosis_functions.scan_diagnosis(CR_data_dict, sub_space_maps_loaded['timeseries'],
                                   self.inputs.common_anat_ref_file, self.inputs.common_mask_file, temporal_info,
                                   spatial_info, plot_seed_frequencies=self.inputs.plot_seed_frequencies,
                                   brainmap_percent_threshold=self.inputs.brainmap_percent_threshold)

        figure_path = os.path.abspath(filename_split)
        spatial_fig_file = figure_path+f'_spatial_diagnosis.{figure_format}'
        fig2.savefig(spatial_fig_file, bbox_inches='tight')

        if len(temporal_fig_list)==1:
            temporal_fig_file = figure_path+f'_temporal_diagnosis.{figure_format}'
            temporal_fig_list[0].savefig(temporal_fig_file, bbox_inches='tight')
        else: # long scans have figures stacked in PDF
            from matplotlib.backends.backend_pdf import PdfPages
            temporal_fig_file = figure_path+f'_temporal_diagnosis.pdf'
            with PdfPages(temporal_fig_file) as pdf:
                for fig in temporal_fig_list:
                    pdf.savefig(fig, bbox_inches='tight')

        setattr(self, 'figure_temporal_diagnosis',
                temporal_fig_file)
        setattr(self, 'figure_spatial_diagnosis',
                spatial_fig_file)
        setattr(self, 'temporal_info', temporal_info)
        setattr(self, 'spatial_info', spatial_info)

        return runtime

    def _list_outputs(self):
        return {'figure_temporal_diagnosis': getattr(self, 'figure_temporal_diagnosis'),
                'figure_spatial_diagnosis': getattr(self, 'figure_spatial_diagnosis'),
                'temporal_info': getattr(self, 'temporal_info'),
                'spatial_info': getattr(self, 'spatial_info'), }


class DatasetDiagnosisInputSpec(BaseInterfaceInputSpec):
    scan_data_list = traits.List(
        exists=True, mandatory=True, desc="A dictionary regrouping the all required accompanying data per scan.")
    seed_prior_maps = traits.List(
        exists=True, desc="A list of expected network map associated to each seed-FC.")
    outlier_threshold = traits.Float(
        desc="The threshold for modified Z-score to classify outliers.")
    network_weighting = traits.Str(
        desc="Whether network maps are absolute or relative.")
    scan_QC_thresholds = traits.Dict(
        desc="Specifications for scan-level QC thresholds.")
    group_avg_prior = traits.Bool(
        desc="Whether to use the group average (median) as a network prior instead of using an external image.")
    figure_format = traits.Str(
        desc="Select file format for figures.")
    extended_QC = traits.Bool(
        desc="Whether to include image intensity and BOLDsd in the group stats.")
    brainmap_percent_threshold = traits.Float(
        desc="Input percentage value for thresholding images.")
    add_smoothing = traits.Bool(
        desc="When to smooth all brain maps prior to QC metric computation with 0.3mm kernel.")
    interpolation = traits.Int(
        desc="Provide an SITK interpolator, used to resample prior_maps into commonspace.")
    rabies_data_type = traits.Int(
        desc="Integer specifying SimpleITK data type.")
    prior_bold_idx = traits.List(
        desc="The index for the ICA components that correspond to bold sources, matching the DR/NPR network ordering.")


class DatasetDiagnosisOutputSpec(TraitedSpec):
    dataset_diagnosis_folder = traits.Str(
        exists=True, desc="Output folder from the dataset diagnosis.")


class DatasetDiagnosis(BaseInterface):
    """
    Conducts a group-level correlation analysis to assess artefact effects.
    Computes the voxelwise cross-subject correlation between each spatial features
    from the previously run scan diagnosis.
    """

    input_spec = DatasetDiagnosisInputSpec
    output_spec = DatasetDiagnosisOutputSpec

    def _run_interface(self, runtime):
        import pathlib
        import matplotlib.pyplot as plt
        from rabies.utils import flatten_list
        from .dataset_QC import generate_dataset_QC,QC_distributions
        from ..utils import load_resample_analysis_maps

        figure_format = self.inputs.figure_format

        out_dir_global = os.path.abspath('dataset_diagnosis_folder/')
        os.makedirs(out_dir_global, exist_ok=True)

        # the output is always the content of this folder
        setattr(self, 'dataset_diagnosis_folder',
                out_dir_global)

        out_dir_parametric = out_dir_global+'/parametric_stats/'
        os.makedirs(out_dir_parametric, exist_ok=True)
        out_dir_non_parametric = out_dir_global+'/non_parametric_stats/'
        os.makedirs(out_dir_non_parametric, exist_ok=True)
        out_dir_dist = out_dir_global+'/sample_distributions/'
        os.makedirs(out_dir_dist, exist_ok=True)

        merged = flatten_list(list(self.inputs.scan_data_list))
        if len(merged) < 3:
            from nipype import logging
            log = logging.getLogger('nipype.workflow')
            log.warning(
                "Cannot run statistics on a sample size smaller than 3, so dataset diagnosis is not run.")
            return runtime

        template_file = merged[0]['anat_ref_file']
        mask_file = merged[0]['mask_file']
        prior_maps_file = merged[0]['prior_maps_file']

        commonspace_maps_loaded = load_resample_analysis_maps(
            mask_file, template_file,
            prior_maps=prior_maps_file,
            interpolation=self.inputs.interpolation, rabies_data_type=self.inputs.rabies_data_type,
            output_prefix=f'dataset_diagnosis_common')
        volume_indices = commonspace_maps_loaded['volume_indices']
        prior_map_vectors = commonspace_maps_loaded['prior_map_vectors']
        if prior_map_vectors is not None:
            prior_bold_idx = [int(i) for i in self.inputs.prior_bold_idx]
            prior_map_vectors = prior_map_vectors[prior_bold_idx]


        scan_name_list=[]
        mean_maps=[]
        std_maps=[]
        CRsd_maps=[]
        FD_DVARS_corr_list=[]
        tdof_list=[]
        mean_FD_list=[]
        total_CRsd_list=[]

        FC_maps_dict={}
        FC_maps_dict['DR']=[]
        FC_maps_dict['NPR']=[]
        FC_maps_dict['SBC']=[]
        
        DR_conf_corr_dict={}
        DR_conf_corr_dict['DR']=[]
        DR_conf_corr_dict['NPR']=[]
        DR_conf_corr_dict['SBC']=[]

        for scan_data in merged:
            temporal_std = scan_data['temporal_std']
            CRsd = scan_data['predicted_std']
            if np.median(CRsd)==0:
                # exclude scans with a majority of 0s; indicates misregistration
                continue
            scan_name = pathlib.Path(scan_data['name_source']).name.rsplit(".nii")[0]
            scan_name_list.append(scan_name)
            mean_maps.append(scan_data['voxelwise_mean'])
            std_maps.append(temporal_std)
            CRsd_maps.append(CRsd)
            total_CRsd_list.append(scan_data['CR_global_std'])
            FD_DVARS_corr_list.append(scan_data['FD_DVARS_corr'])
            tdof_list.append(scan_data['tDOF'])
            mean_FD_list.append(scan_data['FD_trace'].to_numpy().mean())

            if scan_data['DR_bold'] is not None:
                FC_maps_dict['DR'].append(scan_data['DR_bold'])
            if scan_data['NPR_maps'] is not None:
                FC_maps_dict['NPR'].append(scan_data['NPR_maps'])
            if scan_data['seed_map_list'] is not None:
                FC_maps_dict['SBC'].append(scan_data['seed_map_list'])

            # computing the temporal correlation between network and confound timecourses
            DR_confound_time = scan_data['DR_confound_time']
            if DR_confound_time is not None:
                for network_time,key in zip(
                    [scan_data['DR_network_time'],scan_data['NPR_network_time'],scan_data['SBC_network_time']],
                    ['DR','NPR','SBC']):
                    if network_time is not None:
                        # for each network, compute its confound correlation as mean across all DR confound components
                        corr_list = [np.abs(np.corrcoef(network_time[:,[i]].T,DR_confound_time.T)[0,1:]).mean() for i in range(network_time.shape[1])]
                        DR_conf_corr_dict[key].append(corr_list)
 
        # save the list of the scan names that were included in the group statistics
        pd.DataFrame(scan_name_list).to_csv(f'{out_dir_global}/data_diagnosis_scanlist.txt', index=None, header=False)

        from rabies.utils import recover_3D
        std_maps=np.array(std_maps)
        non_zero_voxels = ((std_maps==0).sum(axis=0).astype(bool)==0)
        non_zero_mask = os.path.abspath('non_zero_mask.nii.gz')
        sitk.WriteImage(recover_3D(mask_file, non_zero_voxels.astype(float)), non_zero_mask)

        CRsd_maps=np.array(CRsd_maps)[:,non_zero_voxels]

        corr_variable = []
        variable_name = []
        if self.inputs.extended_QC:
            mean_maps=np.array(mean_maps)[:,non_zero_voxels]
            BOLD_std_maps=np.array(std_maps)[:,non_zero_voxels]
            corr_variable += [mean_maps,BOLD_std_maps]
            variable_name += ['BOLD mean', '$\mathregular{BOLD_{SD}}$']

        corr_variable += [CRsd_maps, np.array(mean_FD_list).reshape(-1,1)]
        variable_name += ['$\mathregular{CR_{SD}}$', 'Mean FD']

        mean_FD_array = np.array(mean_FD_list)
        total_CRsd = np.array(total_CRsd_list)
        FD_DVARS_corr = np.array(FD_DVARS_corr_list)

        # tdof effect; if there's no variability don't compute
        if not np.array(tdof_list).std()==0:
            tdof = np.array(tdof_list).reshape(-1,1)
            corr_variable.append(tdof)
            variable_name.append('tDOF')
            tdof_array = np.array(tdof_list)
        else:
            tdof_array = None

        def change_columns(df):
            columns = list(df.columns)
            i=0
            for column in columns:
                if '$\mathregular{CR_{SD}}$' in column:
                    if 'Overlap:' in column:
                        columns[i] = 'Overlap: Prior - CRsd'
                    if 'Avg.:' in column:
                        columns[i] = 'Avg.: CRsd'
                elif '$\mathregular{BOLD_{SD}}$' in column:
                    if 'Overlap:' in column:
                        columns[i] = 'Overlap: Prior - BOLDsd'
                    if 'Avg.:' in column:
                        columns[i] = 'Avg.: BOLDsd'
                i+=1
            df.columns = columns
            return df
                        

        def generate_dataset_QC_network_i(i,FC_maps,prior_map,non_zero_mask, corr_variable, variable_name, template_file, out_dir_parametric, out_dir_non_parametric,analysis_prefix):

            for non_parametric,out_dir in zip([False, True], [out_dir_parametric, out_dir_non_parametric]):
                with np.errstate(invalid='ignore', divide='ignore'):
                    dataset_stats,fig,fig_unthresholded = generate_dataset_QC(FC_maps, prior_map, non_zero_mask, corr_variable, variable_name, template_file, non_parametric=non_parametric, top_percent=self.inputs.brainmap_percent_threshold, smoothing=self.inputs.add_smoothing)
                df = pd.DataFrame(dataset_stats, index=[1])
                df = change_columns(df)
                df.to_csv(f'{out_dir}/{analysis_prefix}{i}_QC_stats.csv', index=None)
                fig_path = f'{out_dir}/{analysis_prefix}{i}_QC_maps.{figure_format}'
                fig.savefig(fig_path, bbox_inches='tight')
                fig_path = f'{out_dir}/{analysis_prefix}{i}_QC_maps_unthresholded.{figure_format}'
                fig_unthresholded.savefig(fig_path, bbox_inches='tight')

                plt.close(fig)
                plt.close(fig_unthresholded)

        def distribution_network_i(i,prior_map,FC_maps,network_var,DR_conf_corr,FD_DVARS_corr, total_CRsd, mean_FD_array, tdof_array, scan_name_list, outlier_threshold,out_dir_dist,scan_QC_thresholds, analysis_prefix):
            ### PLOT DISTRIBUTIONS FOR OUTLIER DETECTION
            with np.errstate(invalid='ignore', divide='ignore'):
                fig,df,QC_inclusion = QC_distributions(prior_map,FC_maps,network_var,DR_conf_corr,FD_DVARS_corr,total_CRsd, mean_FD_array, tdof_array, scan_name_list, scan_QC_thresholds=scan_QC_thresholds, outlier_threshold=outlier_threshold, top_percent=self.inputs.brainmap_percent_threshold)
            df.to_csv(f'{out_dir_dist}/{analysis_prefix}{i}_outlier_detection.csv', index=None)
            fig_path = f'{out_dir_dist}/{analysis_prefix}{i}_sample_distribution.{figure_format}'
            fig.savefig(fig_path, bbox_inches='tight')
            plt.close(fig)
            return QC_inclusion
        
        def prep_QC_thresholds_i(scan_QC_thresholds, analysis, network_i, num_priors):
            analysis_keys = list(scan_QC_thresholds.keys())
            if not analysis in analysis_keys:
                return {'Dice':None, 'Conf':None, 'Amp':False}
            QC_sub_dict = scan_QC_thresholds[analysis]
            QC_thresholds_i={}
            keys = list(QC_sub_dict.keys())
            for key in ['Dice','Conf']:
                if key in keys:
                    if not type(QC_sub_dict[key]) is list:
                        raise ValueError(f"'{QC_sub_dict[key]}' must be a list.")
                    elif len(QC_sub_dict[key])==0:
                        QC_thresholds_i[key]=None
                    else:
                        if not len(QC_sub_dict[key])==num_priors:
                            raise ValueError(f"The number of Dice thresholds for --scan_QC_thresholds does not match the number of {analysis} networks.")
                        QC_thresholds_i[key]=QC_sub_dict[key][network_i]
                else:
                    QC_thresholds_i[key]=None

            if 'Amp' in keys:
                if not type(QC_sub_dict['Amp']) is bool:
                    raise ValueError(f"'{QC_sub_dict['Amp']}' must be True or False.")
                QC_thresholds_i['Amp']=QC_sub_dict['Amp']
            else:
                QC_thresholds_i['Amp']=False
            return QC_thresholds_i

        scan_QC_thresholds = self.inputs.scan_QC_thresholds


        if len(FC_maps_dict['DR'])>0:
            DR_maps_list=np.array(FC_maps_dict['DR'])

            if self.inputs.group_avg_prior:
                num_priors = DR_maps_list.shape[1]
                prior_maps = np.median(DR_maps_list,axis=0)[:,non_zero_voxels]
            else:
                prior_maps = prior_map_vectors[:,non_zero_voxels]
                num_priors = prior_maps.shape[0]

            for i in range(num_priors):
                if self.inputs.network_weighting=='relative':
                    network_var=None
                else:
                    # network amplitude as L2-norm of a connectivity map
                    network_var = np.sqrt((DR_maps_list[:,i,:] ** 2).sum(axis=1))
                DR_i_scan_QC_thresholds=prep_QC_thresholds_i(scan_QC_thresholds, analysis='DR', network_i=i, num_priors=num_priors)

                FC_maps = DR_maps_list[:,i,non_zero_voxels]
                DR_conf_corr = np.array(DR_conf_corr_dict['DR'])[:,i] if len(DR_conf_corr_dict['DR'])>0 else None
                QC_inclusion = distribution_network_i(i,prior_maps[i,:],FC_maps,network_var,DR_conf_corr,FD_DVARS_corr,total_CRsd, mean_FD_array, tdof_array, scan_name_list, self.inputs.outlier_threshold, out_dir_dist,scan_QC_thresholds=DR_i_scan_QC_thresholds, analysis_prefix='DR')

                # compute group stats only if there is at least 3 scans
                if QC_inclusion.sum()>2:
                    # apply QC inclusion
                    FC_maps_ = FC_maps[QC_inclusion,:]
                    corr_variable_ = [var[QC_inclusion,:] for var in corr_variable]

                    generate_dataset_QC_network_i(i,FC_maps_,prior_maps[i,:],non_zero_mask, corr_variable_, variable_name, template_file, out_dir_parametric, out_dir_non_parametric, analysis_prefix='DR')

        if len(FC_maps_dict['NPR'])>0:
            NPR_maps_list=np.array(FC_maps_dict['NPR'])
            if self.inputs.group_avg_prior:
                num_priors = NPR_maps_list.shape[1]
                prior_maps = np.median(NPR_maps_list,axis=0)[:,non_zero_voxels]
            else:
                prior_maps = prior_map_vectors[:,non_zero_voxels]
                num_priors = prior_maps.shape[0]

            for i in range(num_priors):
                if self.inputs.network_weighting=='relative':
                    network_var=None
                else:
                    # network amplitude as L2-norm of a connectivity map
                    network_var = np.sqrt((NPR_maps_list[:,i,:] ** 2).sum(axis=1))

                NPR_i_scan_QC_thresholds=prep_QC_thresholds_i(scan_QC_thresholds, analysis='NPR', network_i=i, num_priors=num_priors)

                FC_maps = NPR_maps_list[:,i,non_zero_voxels]
                DR_conf_corr = np.array(DR_conf_corr_dict['NPR'])[:,i] if len(DR_conf_corr_dict['NPR'])>0 else None
                QC_inclusion = distribution_network_i(i,prior_maps[i,:],FC_maps,network_var,DR_conf_corr,FD_DVARS_corr,total_CRsd, mean_FD_array, tdof_array, scan_name_list, self.inputs.outlier_threshold, out_dir_dist,scan_QC_thresholds=NPR_i_scan_QC_thresholds, analysis_prefix='NPR')

                # compute group stats only if there is at least 3 scans
                if QC_inclusion.sum()>2:
                    # apply QC inclusion
                    FC_maps_ = FC_maps[QC_inclusion,:]
                    corr_variable_ = [var[QC_inclusion,:] for var in corr_variable]

                    generate_dataset_QC_network_i(i,FC_maps_,prior_maps[i,:],non_zero_mask, corr_variable_, variable_name, template_file, out_dir_parametric, out_dir_non_parametric, analysis_prefix='NPR')

        if len(FC_maps_dict['SBC'])>0:
            seed_maps_list=np.array(FC_maps_dict['SBC'])

            if self.inputs.group_avg_prior:
                num_priors = seed_maps_list.shape[1]
                prior_maps = np.median(seed_maps_list,axis=0)[:,non_zero_voxels]

            # prior maps are provided for seed-FC, tries to run the diagnosis on seeds
            elif len(self.inputs.seed_prior_maps)>0:
                prior_maps=[]
                for prior_map in self.inputs.seed_prior_maps:
                    # resample to match the subject
                    sitk_img = sitk.Resample(sitk.ReadImage(prior_map), sitk.ReadImage(mask_file))
                    prior_maps.append(sitk.GetArrayFromImage(sitk_img)[volume_indices])
                prior_maps = np.array(prior_maps)[:,non_zero_voxels]
                num_priors = prior_maps.shape[0]
            else:
                raise ValueError("Must select either --group_avg_prior or --seed_prior_maps as an option for group statistics with SBC.")
            
            for i in range(num_priors):
                # network amplitude as L2-norm of a connectivity map
                network_var = np.sqrt((seed_maps_list[:,i,:] ** 2).sum(axis=1))

                SBC_i_scan_QC_thresholds=prep_QC_thresholds_i(scan_QC_thresholds, analysis='SBC', network_i=i, num_priors=num_priors)

                FC_maps = seed_maps_list[:,i,non_zero_voxels]
                DR_conf_corr = np.array(DR_conf_corr_dict['SBC'])[:,i] if len(DR_conf_corr_dict['SBC'])>0 else None
                QC_inclusion = distribution_network_i(i,prior_maps[i,:],FC_maps,network_var,DR_conf_corr,FD_DVARS_corr,total_CRsd, mean_FD_array, tdof_array, scan_name_list, self.inputs.outlier_threshold, out_dir_dist,scan_QC_thresholds=SBC_i_scan_QC_thresholds, analysis_prefix='seed_FC')

                # compute group stats only if there is at least 3 scans
                if QC_inclusion.sum()>2:
                    # apply QC inclusion
                    FC_maps_ = FC_maps[QC_inclusion,:]
                    corr_variable_ = [var[QC_inclusion,:] for var in corr_variable]

                    generate_dataset_QC_network_i(i,FC_maps_,prior_maps[i,:],non_zero_mask, corr_variable_, variable_name, template_file, out_dir_parametric, out_dir_non_parametric, analysis_prefix='seed_FC')

        return runtime

    def _list_outputs(self):
        return {
            'dataset_diagnosis_folder': getattr(self, 'dataset_diagnosis_folder'),
            }

