import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def init_bold_confs_wf(opts, aCompCor_method='50%', name="bold_confs_wf"):

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold', 'ref_bold', 'movpar_file', 't1_mask', 't1_labels', 'WM_mask', 'CSF_mask', 'vascular_mask', 'name_source']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['cleaned_bold', 'GSR_cleaned_bold', 'brain_mask', 'WM_mask', 'CSF_mask', 'EPI_labels', 'confounds_csv', 'FD_csv', 'FD_voxelwise', 'pos_voxelwise']),
        name='outputnode')

    WM_mask_to_EPI = pe.Node(MaskEPI(), name='WM_mask_EPI')
    WM_mask_to_EPI.inputs.name_spec = 'WM_mask'

    CSF_mask_to_EPI = pe.Node(MaskEPI(), name='CSF_mask_EPI')
    CSF_mask_to_EPI.inputs.name_spec = 'CSF_mask'

    vascular_mask_to_EPI = pe.Node(MaskEPI(), name='vascular_mask_EPI')
    vascular_mask_to_EPI.inputs.name_spec = 'vascular_mask'

    brain_mask_to_EPI = pe.Node(MaskEPI(), name='Brain_mask_EPI')
    brain_mask_to_EPI.inputs.name_spec = 'brain_mask'

    propagate_labels = pe.Node(MaskEPI(), name='prop_labels_EPI')
    propagate_labels.inputs.name_spec = 'anat_labels'

    estimate_confounds = pe.Node(EstimateConfounds(aCompCor_method=aCompCor_method, rabies_data_type=opts.data_type),
                                 name='estimate_confounds', mem_gb=2.3*opts.scale_min_memory)
    estimate_confounds.plugin_args = {
        'qsub_args': '-pe smp %s' % (str(2*opts.min_proc)), 'overwrite': True}

    workflow = pe.Workflow(name=name)
    workflow.connect([
        (inputnode, WM_mask_to_EPI, [
            ('WM_mask', 'mask'),
            ('ref_bold', 'ref_EPI'),
            ('name_source', 'name_source')]),
        (inputnode, CSF_mask_to_EPI, [
            ('CSF_mask', 'mask'),
            ('ref_bold', 'ref_EPI'),
            ('name_source', 'name_source')]),
        (inputnode, vascular_mask_to_EPI, [
            ('vascular_mask', 'mask'),
            ('ref_bold', 'ref_EPI'),
            ('name_source', 'name_source')]),
        (inputnode, brain_mask_to_EPI, [
            ('t1_mask', 'mask'),
            ('ref_bold', 'ref_EPI'),
            ('name_source', 'name_source')]),
        (inputnode, propagate_labels, [
            ('t1_labels', 'mask'),
            ('ref_bold', 'ref_EPI'),
            ('name_source', 'name_source')]),
        (inputnode, estimate_confounds, [
            ('movpar_file', 'movpar_file'),
            ]),
        (inputnode, estimate_confounds, [
            ('bold', 'bold'),
            ]),
        (WM_mask_to_EPI, estimate_confounds, [
            ('EPI_mask', 'WM_mask')]),
        (WM_mask_to_EPI, outputnode, [
            ('EPI_mask', 'WM_mask')]),
        (CSF_mask_to_EPI, estimate_confounds, [
            ('EPI_mask', 'CSF_mask')]),
        (CSF_mask_to_EPI, outputnode, [
            ('EPI_mask', 'CSF_mask')]),
        (vascular_mask_to_EPI, estimate_confounds, [
            ('EPI_mask', 'vascular_mask')]),
        (brain_mask_to_EPI, estimate_confounds, [
            ('EPI_mask', 'brain_mask')]),
        (brain_mask_to_EPI, outputnode, [
            ('EPI_mask', 'brain_mask')]),
        (propagate_labels, outputnode, [
            ('EPI_mask', 'EPI_labels')]),
        (estimate_confounds, outputnode, [
            ('confounds_csv', 'confounds_csv'),
            ('FD_csv', 'FD_csv'),
            ('FD_voxelwise', 'FD_voxelwise'),
            ('pos_voxelwise', 'pos_voxelwise'),
            ]),
        ])

    return workflow


class EstimateConfoundsInputSpec(BaseInterfaceInputSpec):
    bold = File(exists=True, mandatory=True,
                desc="Preprocessed bold file to clean")
    movpar_file = File(exists=True, mandatory=True,
                       desc="CSV file with the 6 rigid body parameters")
    brain_mask = File(exists=True, mandatory=True,
                      desc="EPI-formated whole brain mask")
    WM_mask = File(exists=True, mandatory=True,
                   desc="EPI-formated white matter mask")
    CSF_mask = File(exists=True, mandatory=True, desc="EPI-formated CSF mask")
    vascular_mask = File(exists=True, mandatory=True,
                         desc="EPI-formated vascular mask")
    aCompCor_method = traits.Str(
        desc="The type of evaluation for the number of aCompCor components: either '50%' or 'first_5'.")
    rabies_data_type = traits.Int(mandatory=True,
        desc="Integer specifying SimpleITK data type.")


class EstimateConfoundsOutputSpec(TraitedSpec):
    confounds_csv = traits.File(desc="CSV file of confounds")
    FD_csv = traits.File(desc="CSV file with global framewise displacement.")
    FD_voxelwise = traits.File(
        desc=".nii file with voxelwise framewise displacement.")
    pos_voxelwise = traits.File(desc=".nii file with voxelwise Positioning.")


class EstimateConfounds(BaseInterface):

    input_spec = EstimateConfoundsInputSpec
    output_spec = EstimateConfoundsOutputSpec

    def _run_interface(self, runtime):
        import numpy as np
        import os
        from rabies.preprocess_pkg.utils import run_command
        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(self.inputs.bold).name.rsplit(".nii")

        # generate a .nii file representing the positioning or framewise displacement for each voxel within the brain_mask
        # first the voxelwise positioning map
        command = 'antsMotionCorrStats -m %s -o %s_pos_file.csv -x %s \
                    -d %s -s %s_pos_voxelwise.nii.gz' % (self.inputs.movpar_file, filename_split[0], self.inputs.brain_mask, self.inputs.bold, filename_split[0])
        rc = run_command(command)
        pos_voxelwise = os.path.abspath(
            "%s_pos_file.nii.gz" % filename_split[0])

        # then the voxelwise framewise displacement map
        command = 'antsMotionCorrStats -m %s -o %s_FD_file.csv -x %s \
                    -d %s -s %s_FD_voxelwise.nii.gz -f 1' % (self.inputs.movpar_file, filename_split[0], self.inputs.brain_mask, self.inputs.bold, filename_split[0])
        rc = run_command(command)

        FD_csv = os.path.abspath("%s_FD_file.csv" % filename_split[0])
        FD_voxelwise = os.path.abspath("%s_FD_file.nii.gz" % filename_split[0])

        confounds = []
        csv_columns = []
        WM_signal = extract_mask_trace(self.inputs.bold, self.inputs.WM_mask)
        confounds.append(WM_signal)
        csv_columns += ['WM_signal']

        CSF_signal = extract_mask_trace(self.inputs.bold, self.inputs.CSF_mask)
        confounds.append(CSF_signal)
        csv_columns += ['CSF_signal']

        vascular_signal = extract_mask_trace(
            self.inputs.bold, self.inputs.vascular_mask)
        confounds.append(vascular_signal)
        csv_columns += ['vascular_signal']

        [aCompCor, num_comp] = compute_aCompCor(
            self.inputs.bold, self.inputs.WM_mask, self.inputs.CSF_mask, method=self.inputs.aCompCor_method, rabies_data_type=self.inputs.rabies_data_type)
        for param in range(aCompCor.shape[1]):
            confounds.append(aCompCor[:, param])
        comp_column = []
        for comp in range(num_comp):
            comp_column.append('aCompCor'+str(comp+1))
        csv_columns += comp_column

        global_signal = extract_mask_trace(
            self.inputs.bold, self.inputs.brain_mask)
        confounds.append(global_signal)
        csv_columns += ['global_signal']
        motion_24 = motion_24_params(self.inputs.movpar_file)
        for param in range(motion_24.shape[1]):
            confounds.append(motion_24[:, param])
        csv_columns += ['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3', 'mov1_der', 'mov2_der', 'mov3_der', 'rot1_der', 'rot2_der', 'rot3_der',
                        'mov1^2', 'mov2^2', 'mov3^2', 'rot1^2', 'rot2^2', 'rot3^2', 'mov1_der^2', 'mov2_der^2', 'mov3_der^2', 'rot1_der^2', 'rot2_der^2', 'rot3_der^2']

        confounds_csv = write_confound_csv(np.transpose(
            np.asarray(confounds)), csv_columns, filename_split[0])

        setattr(self, 'FD_csv', FD_csv)
        setattr(self, 'FD_voxelwise', FD_voxelwise)
        setattr(self, 'pos_voxelwise', pos_voxelwise)
        setattr(self, 'confounds_csv', confounds_csv)
        return runtime

    def _list_outputs(self):
        return {'confounds_csv': getattr(self, 'confounds_csv'),
                'FD_csv': getattr(self, 'FD_csv'),
                'pos_voxelwise': getattr(self, 'pos_voxelwise'),
                'FD_voxelwise': getattr(self, 'FD_voxelwise')}


def write_confound_csv(confound_array, column_names, filename_template):
    import pandas as pd
    df = pd.DataFrame(confound_array)
    df.columns = column_names
    csv_path = os.path.abspath("%s_confounds.csv" % filename_template)
    df.to_csv(csv_path)
    return csv_path


def compute_aCompCor(bold, WM_mask, CSF_mask, method='50%', rabies_data_type=8):
    '''
    Compute the anatomical comp corr through PCA over a defined ROI (mask) within
    the EPI, and retain either the first 5 components' time series or up to 50% of
    the variance explained as in Muschelli et al. 2014.
    '''
    import SimpleITK as sitk
    from sklearn.decomposition import PCA

    WM_data = sitk.GetArrayFromImage(sitk.ReadImage(
        WM_mask, rabies_data_type))
    CSF_data = sitk.GetArrayFromImage(sitk.ReadImage(
        CSF_mask, rabies_data_type))
    combined = (WM_data+CSF_data) > 0
    from .utils import copyInfo_3DImage
    noise_mask = copyInfo_3DImage(sitk.GetImageFromArray(combined.astype(
        int), isVector=False), sitk.ReadImage(WM_mask, sitk.sitkInt16))
    sitk.WriteImage(noise_mask, 'noise_mask.nii.gz')

    from nilearn.input_data import NiftiMasker
    # detrend and standardize the voxel time series before PCA
    masker = NiftiMasker(mask_img='noise_mask.nii.gz',
                         standardize=True, detrend=True)
    mask_timeseries = masker.fit_transform(
        bold)  # shape n_timepoints x n_voxels

    if method == '50%':
        pca = PCA()
        pca.fit(mask_timeseries)
        explained_variance = pca.explained_variance_ratio_
        cum_var = 0
        num_comp = 0
        # evaluate the # of components to explain 50% of the variance
        while(cum_var <= 0.5):
            cum_var += explained_variance[num_comp]
            num_comp += 1
    elif method == 'first_5':
        num_comp = 5

    pca = PCA(n_components=num_comp)
    comp_timeseries = pca.fit_transform(mask_timeseries)
    import logging
    log = logging.getLogger('root')
    log.debug("Extracting "+str(num_comp)+" components for aCompCorr.")
    return comp_timeseries, num_comp


def motion_24_params(movpar_csv):
    '''
    motioncorr_24params: 6 head motion parameters, their temporal derivative, and the 12 corresponding squared items (Friston et al. 1996, Magn. Reson. Med.)
    '''
    import numpy as np
    rigid_params = extract_rigid_movpar(movpar_csv)
    movpar = np.zeros([np.size(rigid_params, 0), 24])
    movpar[:, :6] = rigid_params
    for i in range(6):
        # Compute temporal derivative as difference between two neighboring points
        movpar[0, 6+i] = 0
        movpar[1:, 6+i] = movpar[1:, i]-movpar[:-1, i]
        # add the squared coefficients
        movpar[:, 12+i] = movpar[:, i]**2
        movpar[:, 18+i] = movpar[:, 6+i]**2
    return movpar


def extract_rigid_movpar(movpar_csv):
    import numpy as np
    import csv
    temp = []
    with open(movpar_csv) as csvfile:
        motcorr = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in motcorr:
            temp.append(row)
    movpar = np.zeros([(len(temp)-1), 6])
    j = 0
    for row in temp[1:]:
        for i in range(2, len(row)):
            movpar[j, i-2] = float(row[i])
        j = j+1
    return movpar


def extract_mask_trace(bold, mask):
    import numpy as np
    import nilearn.masking
    mask_signal = nilearn.masking.apply_mask(bold, mask)
    mean_trace = np.mean(mask_signal, 1)
    return mean_trace


def extract_labels(atlas):
    import nilearn.regions
    nilearn.regions.connected_label_regions(atlas)


class MaskEPIInputSpec(BaseInterfaceInputSpec):
    mask = File(exists=True, mandatory=True,
                desc="Mask to transfer to EPI space.")
    ref_EPI = File(exists=True, mandatory=True,
                   desc="Motion-realigned and SDC-corrected reference 3D EPI.")
    name_spec = traits.Str(desc="Specify the name of the mask.")
    name_source = File(exists=True, mandatory=True,
                       desc='Reference BOLD file for naming the output.')


class MaskEPIOutputSpec(TraitedSpec):
    EPI_mask = traits.File(desc="The generated EPI mask.")


class MaskEPI(BaseInterface):

    input_spec = MaskEPIInputSpec
    output_spec = MaskEPIOutputSpec

    def _run_interface(self, runtime):
        import os
        import SimpleITK as sitk

        import pathlib  # Better path manipulation
        filename_split = pathlib.Path(
            self.inputs.name_source).name.rsplit(".nii")

        if self.inputs.name_spec is None:
            new_mask_path = os.path.abspath(
                '%s_EPI_mask.nii.gz' % (filename_split[0],))
        else:
            new_mask_path = os.path.abspath('%s_%s.nii.gz' % (
                filename_split[0], self.inputs.name_spec,))

        command = 'antsApplyTransforms -i ' + self.inputs.mask + ' -r ' + \
            self.inputs.ref_EPI + ' -o ' + new_mask_path + ' -n GenericLabel'
        from rabies.preprocess_pkg.utils import run_command
        rc = run_command(command)

        sitk.WriteImage(sitk.ReadImage(
            new_mask_path, sitk.sitkInt16), new_mask_path)

        setattr(self, 'EPI_mask', new_mask_path)
        return runtime

    def _list_outputs(self):
        return {'EPI_mask': getattr(self, 'EPI_mask')}
