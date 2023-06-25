import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
from rabies.analysis_pkg.diagnosis_pkg import diagnosis_functions

from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)

class PrepMasksInputSpec(BaseInterfaceInputSpec):
    mask_dict_list = traits.List(
        exists=True, mandatory=True, desc="Brain mask.")
    prior_maps = File(exists=True, mandatory=True,
                      desc="MELODIC ICA components to use.")
    DSURQE_regions = traits.Bool(
        desc="Whether to use the regional masks generated from the DSURQE atlas for the grayplots outputs. Requires using the DSURQE template for preprocessing.")


class PrepMasksOutputSpec(TraitedSpec):
    mask_file_dict = traits.Dict(
        desc="A dictionary regrouping the all required accompanying files.")


class PrepMasks(BaseInterface):
    """

    """

    input_spec = PrepMasksInputSpec
    output_spec = PrepMasksOutputSpec

    def _run_interface(self, runtime):
        from rabies.utils import flatten_list,resample_image_spacing
        merged = flatten_list(list(self.inputs.mask_dict_list))
        mask_dict = merged[0]  # all mask files are assumed to be identical
        brain_mask_file = mask_dict['mask_file']
        WM_mask_file = mask_dict['WM_mask_file']
        CSF_mask_file = mask_dict['CSF_mask_file']

        # resample the template to the EPI dimensions
        resampled = resample_image_spacing(sitk.ReadImage(mask_dict['preprocess_anat_template']), sitk.ReadImage(
            brain_mask_file).GetSpacing())
        template_file = os.path.abspath('display_template.nii.gz')
        sitk.WriteImage(resampled, template_file)

        if self.inputs.DSURQE_regions:
            if 'XDG_DATA_HOME' in os.environ.keys():
                rabies_path = os.environ['XDG_DATA_HOME']+'/rabies'
            else:
                rabies_path = os.environ['HOME']+'/.local/share/rabies'
            right_hem_mask_file = diagnosis_functions.resample_mask(rabies_path+'/DSURQE_40micron_right_hem_mask.nii.gz',
                                                brain_mask_file)
            left_hem_mask_file = diagnosis_functions.resample_mask(rabies_path+'/DSURQE_40micron_left_hem_mask.nii.gz',
                                               brain_mask_file)
        else:
            right_hem_mask_file = ''
            left_hem_mask_file = ''

        from rabies.analysis_pkg.analysis_functions import resample_IC_file
        prior_maps = resample_IC_file(self.inputs.prior_maps, brain_mask_file)

        edge_mask_file = os.path.abspath('edge_mask.nii.gz')
        diagnosis_functions.compute_edge_mask(brain_mask_file, edge_mask_file, num_edge_voxels=1)
        mask_file_dict = {'template_file': template_file, 'brain_mask': brain_mask_file, 'WM_mask': WM_mask_file, 'CSF_mask': CSF_mask_file,
                          'edge_mask': edge_mask_file, 'right_hem_mask': right_hem_mask_file, 'left_hem_mask': left_hem_mask_file, 'prior_maps': prior_maps}

        setattr(self, 'mask_file_dict', mask_file_dict)
        return runtime

    def _list_outputs(self):
        return {'mask_file_dict': getattr(self, 'mask_file_dict')}


class ScanDiagnosisInputSpec(BaseInterfaceInputSpec):
    dict_file = File(exists=True, mandatory=True, desc="Dictionary with prepared analysis data.")
    analysis_dict = traits.Dict(
        desc="A dictionary regrouping relevant outputs from analysis.")
    prior_bold_idx = traits.List(
        desc="The index for the ICA components that correspond to bold sources.")
    prior_confound_idx = traits.List(
        desc="The index for the ICA components that correspond to confounding sources.")
    DSURQE_regions = traits.Bool(
        desc="Whether to use the regional masks generated from the DSURQE atlas for the grayplots outputs. Requires using the DSURQE template for preprocessing.")
    figure_format = traits.Str(
        desc="Select file format for figures.")


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
        import pickle
        with open(self.inputs.dict_file, 'rb') as handle:
            data_dict = pickle.load(handle)

        figure_format = self.inputs.figure_format

        # convert to an integer list
        prior_bold_idx = [int(i) for i in self.inputs.prior_bold_idx]
        prior_confound_idx = [int(i) for i in self.inputs.prior_confound_idx]

        temporal_info, spatial_info = diagnosis_functions.process_data(
            data_dict, self.inputs.analysis_dict, prior_bold_idx, prior_confound_idx)

        fig, fig2 = diagnosis_functions.scan_diagnosis(data_dict, temporal_info,
                                   spatial_info, regional_grayplot=self.inputs.DSURQE_regions)

        import pathlib
        filename_template = pathlib.Path(data_dict['name_source']).name.rsplit(".nii")[0]
        figure_path = os.path.abspath(filename_template)
        fig.savefig(figure_path+f'_temporal_diagnosis.{figure_format}', bbox_inches='tight')
        fig2.savefig(figure_path+f'_spatial_diagnosis.{figure_format}', bbox_inches='tight')

        setattr(self, 'figure_temporal_diagnosis',
                figure_path+f'_temporal_diagnosis.{figure_format}')
        setattr(self, 'figure_spatial_diagnosis',
                figure_path+f'_spatial_diagnosis.{figure_format}')
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
    figure_format = traits.Str(
        desc="Select file format for figures.")
    extended_QC = traits.Bool(
        desc="Whether to include image intensity and BOLDsd in the group stats.")


class DatasetDiagnosisOutputSpec(TraitedSpec):
    analysis_QC = traits.Str(
        exists=True, desc="Output figure from the analysis QC.")


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
        from .analysis_QC import analysis_QC,QC_distributions

        figure_format = self.inputs.figure_format

        out_dir_global = os.path.abspath('analysis_QC/')
        os.makedirs(out_dir_global, exist_ok=True)
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
            setattr(self, 'analysis_QC',
                    out_dir_global)
            return runtime

        template_file = merged[0]['template_file']
        mask_file = merged[0]['mask_file']
        brain_mask = sitk.GetArrayFromImage(sitk.ReadImage(mask_file))
        volume_indices = brain_mask.astype(bool)

        scan_name_list=[]
        mean_maps=[]
        std_maps=[]
        CRsd_maps=[]
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
            tdof_list.append(scan_data['tDOF'])
            mean_FD_list.append(scan_data['FD_trace'].to_numpy().mean())

            FC_maps_dict['DR'].append(scan_data['DR_BOLD'])
            FC_maps_dict['NPR'].append(scan_data['NPR_maps'])
            FC_maps_dict['SBC'].append(scan_data['seed_map_list'])

            # computing the temporal correlation between network and confound timecourses
            DR_confound_time = scan_data['DR_confound_time']
            for network_time,key in zip(
                [scan_data['DR_network_time'],scan_data['NPR_network_time'],scan_data['SBC_network_time']],
                ['DR','NPR','SBC']):
                if len(network_time)>0:
                    # for each network, compute its confound correlation as mean across all DR confound components
                    corr_list = [np.abs(np.corrcoef(network_time[:,[i]].T,DR_confound_time.T)[0,1:]).mean() for i in range(network_time.shape[1])]
                    DR_conf_corr_dict[key].append(corr_list)
 
        # save the list of the scan names that were included in the group statistics
        pd.DataFrame(scan_name_list).to_csv(f'{out_dir_global}/analysis_QC_scanlist.txt', index=None, header=False)

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
                        

        def analysis_QC_network_i(i,FC_maps,prior_map,non_zero_mask, corr_variable, variable_name, template_file, out_dir_parametric, out_dir_non_parametric,analysis_prefix):

            for non_parametric,out_dir in zip([False, True], [out_dir_parametric, out_dir_non_parametric]):
                dataset_stats,fig,fig_unthresholded = analysis_QC(FC_maps, prior_map, non_zero_mask, corr_variable, variable_name, template_file, non_parametric=non_parametric)
                df = pd.DataFrame(dataset_stats, index=[1])
                df = change_columns(df)
                df.to_csv(f'{out_dir}/{analysis_prefix}{i}_QC_stats.csv', index=None)
                fig_path = f'{out_dir}/{analysis_prefix}{i}_QC_maps.{figure_format}'
                fig.savefig(fig_path, bbox_inches='tight')
                fig_path = f'{out_dir}/{analysis_prefix}{i}_QC_maps_unthresholded.{figure_format}'
                fig_unthresholded.savefig(fig_path, bbox_inches='tight')

                plt.close(fig)
                plt.close(fig_unthresholded)

        def distribution_network_i(i,prior_map,FC_maps,network_var,DR_conf_corr,total_CRsd, mean_FD_array, tdof_array, scan_name_list, outlier_threshold,out_dir_dist,scan_QC_thresholds, analysis_prefix):
            ### PLOT DISTRIBUTIONS FOR OUTLIER DETECTION
            fig,df,QC_inclusion = QC_distributions(prior_map,FC_maps,network_var,DR_conf_corr,total_CRsd, mean_FD_array, tdof_array, scan_name_list, scan_QC_thresholds=scan_QC_thresholds, outlier_threshold=outlier_threshold)
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

        prior_maps = scan_data['prior_maps'][:,non_zero_voxels]
        num_priors = prior_maps.shape[0]

        DR_maps_list=np.array(FC_maps_dict['DR'])
        for i in range(num_priors):
            if self.inputs.network_weighting=='relative':
                network_var=None
            else:
                # we don't apply the non_zero_voxels mask as it changes the original variance estimate
                network_var = np.sqrt((DR_maps_list[:,i,:] ** 2).sum(axis=1)) # the component variance/scaling is taken from the spatial maps
            DR_i_scan_QC_thresholds=prep_QC_thresholds_i(scan_QC_thresholds, analysis='DR', network_i=i, num_priors=num_priors)

            FC_maps = DR_maps_list[:,i,non_zero_voxels]
            QC_inclusion = distribution_network_i(i,prior_maps[i,:],FC_maps,network_var,np.array(DR_conf_corr_dict['DR'])[:,i],total_CRsd, mean_FD_array, tdof_array, scan_name_list, self.inputs.outlier_threshold, out_dir_dist,scan_QC_thresholds=DR_i_scan_QC_thresholds, analysis_prefix='DR')

            # compute group stats only if there is at least 3 scans
            if QC_inclusion.sum()>2:
                # apply QC inclusion
                FC_maps_ = FC_maps[QC_inclusion,:]
                corr_variable_ = [var[QC_inclusion,:] for var in corr_variable]

                analysis_QC_network_i(i,FC_maps_,prior_maps[i,:],non_zero_mask, corr_variable_, variable_name, template_file, out_dir_parametric, out_dir_non_parametric, analysis_prefix='DR')


        NPR_maps_list=np.array(FC_maps_dict['NPR'])
        if NPR_maps_list.shape[1]>0:
            for i in range(num_priors):
                if self.inputs.network_weighting=='relative':
                    network_var=None
                else:
                    # we don't apply the non_zero_voxels mask as it changes the original variance estimate
                    network_var = np.sqrt((NPR_maps_list[:,i,:] ** 2).sum(axis=1)) # the component variance/scaling is taken from the spatial maps

                NPR_i_scan_QC_thresholds=prep_QC_thresholds_i(scan_QC_thresholds, analysis='NPR', network_i=i, num_priors=num_priors)

                FC_maps = NPR_maps_list[:,i,non_zero_voxels]
                QC_inclusion = distribution_network_i(i,prior_maps[i,:],FC_maps,network_var,np.array(DR_conf_corr_dict['NPR'])[:,i],total_CRsd, mean_FD_array, tdof_array, scan_name_list, self.inputs.outlier_threshold, out_dir_dist,scan_QC_thresholds=NPR_i_scan_QC_thresholds, analysis_prefix='NPR')

                # compute group stats only if there is at least 3 scans
                if QC_inclusion.sum()>2:
                    # apply QC inclusion
                    FC_maps_ = FC_maps[QC_inclusion,:]
                    corr_variable_ = [var[QC_inclusion,:] for var in corr_variable]

                    analysis_QC_network_i(i,FC_maps_,prior_maps[i,:],non_zero_mask, corr_variable_, variable_name, template_file, out_dir_parametric, out_dir_non_parametric, analysis_prefix='NPR')

        # prior maps are provided for seed-FC, tries to run the diagnosis on seeds
        if len(self.inputs.seed_prior_maps)>0:
            prior_maps=[]
            for prior_map in self.inputs.seed_prior_maps:
                # resample to match the subject
                sitk_img = sitk.Resample(sitk.ReadImage(prior_map), sitk.ReadImage(mask_file))
                prior_maps.append(sitk.GetArrayFromImage(sitk_img)[volume_indices])

            prior_maps = np.array(prior_maps)[:,non_zero_voxels]
            num_priors = prior_maps.shape[0]
            seed_maps_list=np.array(FC_maps_dict['SBC'])
            for i in range(num_priors):
                network_var=None

                SBC_i_scan_QC_thresholds=prep_QC_thresholds_i(scan_QC_thresholds, analysis='SBC', network_i=i, num_priors=num_priors)

                FC_maps = seed_maps_list[:,i,non_zero_voxels]
                QC_inclusion = distribution_network_i(i,prior_maps[i,:],FC_maps,network_var,np.array(DR_conf_corr_dict['SBC'])[:,i],total_CRsd, mean_FD_array, tdof_array, scan_name_list, self.inputs.outlier_threshold, out_dir_dist,scan_QC_thresholds=SBC_i_scan_QC_thresholds, analysis_prefix='seed_FC')

                # compute group stats only if there is at least 3 scans
                if QC_inclusion.sum()>2:
                    # apply QC inclusion
                    FC_maps_ = FC_maps[QC_inclusion,:]
                    corr_variable_ = [var[QC_inclusion,:] for var in corr_variable]

                    analysis_QC_network_i(i,FC_maps_,prior_maps[i,:],non_zero_mask, corr_variable_, variable_name, template_file, out_dir_parametric, out_dir_non_parametric, analysis_prefix='seed_FC')

        setattr(self, 'analysis_QC',
                out_dir_global)
        return runtime

    def _list_outputs(self):
        return {
            'analysis_QC': getattr(self, 'analysis_QC'),
            }

