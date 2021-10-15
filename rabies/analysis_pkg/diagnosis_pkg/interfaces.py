import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from rabies.analysis_pkg import analysis_functions
import SimpleITK as sitk
from rabies.analysis_pkg.diagnosis_pkg import diagnosis_functions
from rabies.analysis_pkg.analysis_math import elementwise_corrcoef

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
        from rabies.preprocess_pkg.utils import flatten_list
        merged = flatten_list(list(self.inputs.mask_dict_list))
        mask_dict = merged[0]  # all mask files are assumed to be identical
        brain_mask_file = mask_dict['mask_file']
        WM_mask_file = mask_dict['WM_mask_file']
        CSF_mask_file = mask_dict['CSF_mask_file']

        # resample the template to the EPI dimensions
        from rabies.preprocess_pkg.utils import resample_image_spacing
        resampled = resample_image_spacing(sitk.ReadImage(mask_dict['preprocess_anat_template']), sitk.ReadImage(
            brain_mask_file).GetSpacing(), resampling_interpolation='BSpline')
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

        prior_maps = diagnosis_functions.resample_IC_file(self.inputs.prior_maps, brain_mask_file)

        edge_mask_file = os.path.abspath('edge_mask.nii.gz')
        diagnosis_functions.compute_edge_mask(brain_mask_file, edge_mask_file, num_edge_voxels=1)
        mask_file_dict = {'template_file': template_file, 'brain_mask': brain_mask_file, 'WM_mask': WM_mask_file, 'CSF_mask': CSF_mask_file,
                          'edge_mask': edge_mask_file, 'right_hem_mask': right_hem_mask_file, 'left_hem_mask': left_hem_mask_file, 'prior_maps': prior_maps}

        setattr(self, 'mask_file_dict', mask_file_dict)
        return runtime

    def _list_outputs(self):
        return {'mask_file_dict': getattr(self, 'mask_file_dict')}


class ScanDiagnosisInputSpec(BaseInterfaceInputSpec):
    file_dict = traits.Dict(
        desc="A dictionary regrouping the all required accompanying files.")
    mask_file_dict = traits.Dict(
        desc="A dictionary regrouping the all required accompanying files.")
    prior_bold_idx = traits.List(
        desc="The index for the ICA components that correspond to bold sources.")
    prior_confound_idx = traits.List(
        desc="The index for the ICA components that correspond to confounding sources.")
    dual_ICA = traits.Int(
        desc="number of components to compute from dual ICA.")
    DSURQE_regions = traits.Bool(
        desc="Whether to use the regional masks generated from the DSURQE atlas for the grayplots outputs. Requires using the DSURQE template for preprocessing.")


class ScanDiagnosisOutputSpec(TraitedSpec):
    figure_temporal_diagnosis = File(
        exists=True, desc="Output figure from the scan diagnosis")
    figure_spatial_diagnosis = File(
        exists=True, desc="Output figure from the scan diagnosis")
    temporal_info = traits.Dict(
        desc="A dictionary regrouping the temporal features.")
    spatial_info = traits.Dict(
        desc="A dictionary regrouping the spatial features.")
    VE_file = File(exists=True, mandatory=True,
                   desc="Output the VE file for the datasink.")


class ScanDiagnosis(BaseInterface):
    """
    Extracts several spatial and temporal features on the target scan.
    Spatial features include tSTD, CR-R^2 (variance explained from confound regression),
    correlation maps with global signal/DVARS/FD, and network maps from specified
    BOLD priors at the indices of prior_bold_idx.
    Temporal features include grayplot, 6 motion parameters, framewise displacement,
    DVARS, WM/CSV/edge mask timecourses, CR-R^2, and the average amplitude of BOLD and
    confound components seperately.
    """

    input_spec = ScanDiagnosisInputSpec
    output_spec = ScanDiagnosisOutputSpec

    def _run_interface(self, runtime):
        # convert to an integer list
        bold_file = self.inputs.file_dict['bold_file']
        CR_data_dict = self.inputs.file_dict['CR_data_dict']
        VE_file = self.inputs.file_dict['VE_file']
        prior_bold_idx = [int(i) for i in self.inputs.prior_bold_idx]
        prior_confound_idx = [int(i) for i in self.inputs.prior_confound_idx]

        temporal_info, spatial_info = diagnosis_functions.process_data(
            bold_file, CR_data_dict, VE_file, self.inputs.mask_file_dict, prior_bold_idx, prior_confound_idx, dual_ICA=self.inputs.dual_ICA)

        fig, fig2 = diagnosis_functions.scan_diagnosis(bold_file, self.inputs.mask_file_dict, temporal_info,
                                   spatial_info, CR_data_dict, regional_grayplot=self.inputs.DSURQE_regions)

        import pathlib
        filename_template = pathlib.Path(bold_file).name.rsplit(".nii")[0]
        figure_path = os.path.abspath(filename_template)
        fig.savefig(figure_path+'_temporal_diagnosis.png', bbox_inches='tight')
        fig2.savefig(figure_path+'_spatial_diagnosis.png', bbox_inches='tight')

        setattr(self, 'figure_temporal_diagnosis',
                figure_path+'_temporal_diagnosis.png')
        setattr(self, 'figure_spatial_diagnosis',
                figure_path+'_spatial_diagnosis.png')
        setattr(self, 'temporal_info', temporal_info)
        setattr(self, 'spatial_info', spatial_info)
        setattr(self, 'VE_file', VE_file)

        return runtime

    def _list_outputs(self):
        return {'figure_temporal_diagnosis': getattr(self, 'figure_temporal_diagnosis'),
                'figure_spatial_diagnosis': getattr(self, 'figure_spatial_diagnosis'),
                'temporal_info': getattr(self, 'temporal_info'),
                'spatial_info': getattr(self, 'spatial_info'),
                'VE_file': getattr(self, 'VE_file'), }


class DatasetDiagnosisInputSpec(BaseInterfaceInputSpec):
    spatial_info_list = traits.List(
        exists=True, mandatory=True, desc="A dictionary regrouping the spatial features.")
    mask_file_dict = traits.Dict(
        exists=True, mandatory=True, desc="A dictionary regrouping the all required accompanying files.")


class DatasetDiagnosisOutputSpec(TraitedSpec):
    figure_dataset_diagnosis = File(
        exists=True, desc="Output figure from the dataset diagnosis")


class DatasetDiagnosis(BaseInterface):
    """
    Conducts a group-level correlation analysis to assess artefact effects.
    Computes the voxelwise cross-subject correlation between each spatial features
    from the previously run scan diagnosis.
    """

    input_spec = DatasetDiagnosisInputSpec
    output_spec = DatasetDiagnosisOutputSpec

    def _run_interface(self, runtime):
        from rabies.preprocess_pkg.utils import flatten_list
        import nilearn.plotting
        merged = flatten_list(list(self.inputs.spatial_info_list))
        if len(merged) < 3:
            import logging
            log = logging.getLogger('root')
            log.warning(
                "Cannot run statistics on a sample size smaller than 3, so an empty figure is generated.")
            fig, axes = plt.subplots()
            fig.savefig(os.path.abspath(
                'empty_dataset_diagnosis.png'), bbox_inches='tight')

            setattr(self, 'figure_dataset_diagnosis',
                    os.path.abspath('empty_dataset_diagnosis.png'))
            return runtime

        dict_keys = ['temporal_std', 'VE_spatial', 'GS_corr',
                     'DVARS_corr', 'FD_corr', 'DR_BOLD', 'dual_ICA_maps']

        voxelwise_list = []
        for spatial_info in merged:
            sub_list = [spatial_info[key] for key in dict_keys]
            voxelwise_sub = np.array(sub_list[:5])
            if len(sub_list[6]) > 0:
                voxelwise_sub = np.concatenate(
                    (voxelwise_sub, np.array(sub_list[5]), np.array(sub_list[6])), axis=0)
            else:
                voxelwise_sub = np.concatenate(
                    (voxelwise_sub, np.array(sub_list[5])), axis=0)
            voxelwise_list.append(voxelwise_sub)
            num_DR_maps = len(sub_list[5])
            num_prior_maps = len(sub_list[6])
        voxelwise_array = np.array(voxelwise_list)

        label_name = ['temporal_std', 'VE_spatial',
                      'GS_corr', 'DVARS_corr', 'FD_corr']
        label_name += [f'BOLD Dual Regression map {i}' for i in range(num_DR_maps)]
        label_name += [f'BOLD Dual ICA map {i}' for i in range(num_prior_maps)]

        template_file = self.inputs.mask_file_dict['template_file']
        mask_file = self.inputs.mask_file_dict['brain_mask']
        from rabies.preprocess_pkg.preprocess_visual_QC import plot_3d, otsu_scaling
        scaled = otsu_scaling(template_file)

        ncols = 5
        fig, axes = plt.subplots(nrows=voxelwise_array.shape[1], ncols=ncols, figsize=(
            12*ncols, 2*voxelwise_array.shape[1]))
        for i, x_label in zip(range(voxelwise_array.shape[1]), label_name):
            for j, y_label in zip(range(ncols), label_name[:ncols]):
                ax = axes[i, j]
                if i <= j:
                    ax.axis('off')
                    continue

                X = voxelwise_array[:, i, :]
                Y = voxelwise_array[:, j, :]
                corr = elementwise_corrcoef(X, Y)

                plot_3d([ax], scaled, fig, vmin=0, vmax=1, cmap='gray',
                        alpha=1, cbar=False, num_slices=6, planes=('coronal'))
                analysis_functions.recover_3D(
                    mask_file, corr).to_filename('temp_img.nii.gz')
                sitk_img = sitk.ReadImage('temp_img.nii.gz')
                plot_3d([ax], sitk_img, fig, vmin=-0.7, vmax=0.7, cmap='cold_hot',
                        alpha=1, cbar=True, threshold=0.1, num_slices=6, planes=('coronal'))
                ax.set_title(f'Cross-correlation for {x_label} and {y_label}', fontsize=15, color='white')
        fig.savefig(os.path.abspath('dataset_diagnosis.png'),
                    bbox_inches='tight')

        setattr(self, 'figure_dataset_diagnosis',
                os.path.abspath('dataset_diagnosis.png'))
        return runtime

    def _list_outputs(self):
        return {'figure_dataset_diagnosis': getattr(self, 'figure_dataset_diagnosis')}

