import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)
from nipype import Function
from .utils import Skullstrip

def init_bold_confs_wf(TR, motioncorr_24params=False, apply_GSR=False, name="bold_confs_wf"):
    '''
motioncorr_24params: Regression of 6 head motion parameters and autoregressive
                        models of motion: 6 head motion parameters, 6 head motion parameters one time
                        point before, and the 12 corresponding squared items (Friston et al., 1996)
    '''
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold', 'ref_bold', 'movpar_file', 't1_mask', 't1_labels', 'WM_mask', 'CSF_mask']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['cleaned_bold', 'GSR_cleaned_bold', 'EPI_labels', 'confounds_csv']),
        name='outputnode')

    WM_mask_to_EPI=pe.Node(MaskEPI(use_transforms=False), name='WM_mask_EPI')
    CSF_mask_to_EPI=pe.Node(MaskEPI(use_transforms=False), name='CSF_mask_EPI')
    brain_mask_to_EPI=pe.Node(MaskEPI(use_transforms=False), name='Brain_mask_EPI')
    propagate_labels=pe.Node(MaskEPI(use_transforms=False), name='prop_labels_EPI')
    skullstrip=pe.Node(Skullstrip(), name='skullstrip')

    confound_regression=pe.Node(ConfoundRegression(apply_GSR=False, motioncorr_24params=motioncorr_24params, TR=TR), name='confound_regression')

    if apply_GSR:
        GSR_confound_regression=pe.Node(ConfoundRegression(apply_GSR=True, motioncorr_24params=motioncorr_24params, TR=TR), name='GSR_confound_regression')

    workflow = pe.Workflow(name=name)
    workflow.connect([
        (inputnode, WM_mask_to_EPI, [
            ('WM_mask', 'mask'),
            ('ref_bold', 'ref_EPI')]),
        (inputnode, CSF_mask_to_EPI, [
            ('CSF_mask', 'mask'),
            ('ref_bold', 'ref_EPI')]),
        (inputnode, brain_mask_to_EPI, [
            ('t1_mask', 'mask'),
            ('ref_bold', 'ref_EPI')]),
        (inputnode, propagate_labels, [
            ('t1_labels', 'mask'),
            ('ref_bold', 'ref_EPI')]),
        (inputnode, confound_regression, [
            ('movpar_file', 'movpar_file'),
            ]),
        (inputnode, skullstrip, [
            ('bold', 'in_file')]),
        (brain_mask_to_EPI, skullstrip, [
            ('EPI_mask', 'brain_mask')]),
        (skullstrip, confound_regression, [
            ('skullstrip_brain', 'bold'),
            ]),
        (WM_mask_to_EPI, confound_regression, [
            ('EPI_mask', 'WM_mask')]),
        (CSF_mask_to_EPI, confound_regression, [
            ('EPI_mask', 'CSF_mask')]),
        (brain_mask_to_EPI, confound_regression, [
            ('EPI_mask', 'brain_mask')]),
        (propagate_labels, outputnode, [
            ('EPI_mask', 'EPI_labels')]),
        (confound_regression, outputnode, [
            ('cleaned_bold', 'cleaned_bold'),
            ('confounds_csv', 'confounds_csv'),
            ]),
        ])

    if apply_GSR:
        workflow.connect([
            (inputnode, GSR_confound_regression, [
                ('movpar_file', 'movpar_file'),
                ]),
            (skullstrip, GSR_confound_regression, [
                ('skullstrip_brain', 'bold'),
                ]),
            (WM_mask_to_EPI, GSR_confound_regression, [
                ('EPI_mask', 'WM_mask')]),
            (CSF_mask_to_EPI, GSR_confound_regression, [
                ('EPI_mask', 'CSF_mask')]),
            (brain_mask_to_EPI, GSR_confound_regression, [
                ('EPI_mask', 'brain_mask')]),
            (GSR_confound_regression, outputnode, [
                ('cleaned_bold', 'GSR_cleaned_bold'),
                ]),
            ])

    return workflow

class ConfoundRegressionInputSpec(BaseInterfaceInputSpec):
    bold = File(exists=True, mandatory=True, desc="Preprocessed bold file to clean")
    movpar_file = File(exists=True, mandatory=True, desc="CSV file with the 6 rigid body parameters")
    brain_mask = File(exists=True, mandatory=True, desc="EPI-formated whole brain mask")
    WM_mask = File(exists=True, mandatory=True, desc="EPI-formated white matter mask")
    CSF_mask = File(exists=True, mandatory=True, desc="EPI-formated CSF mask")
    apply_GSR = traits.Bool(mandatory=True, desc="Use global signal regression or not")
    motioncorr_24params = traits.Bool(mandatory=True, desc="Apply 24 parameters motion correction or 6 rigid body correction only")
    TR = traits.Float(mandatory=True, desc="Repetition time.")

class ConfoundRegressionOutputSpec(TraitedSpec):
    cleaned_bold = traits.File(desc="The cleaned bold")
    confounds_csv = traits.File(desc="CSV file of confounds")

class ConfoundRegression(BaseInterface):

    input_spec = ConfoundRegressionInputSpec
    output_spec = ConfoundRegressionOutputSpec

    def _run_interface(self, runtime):
        import numpy as np
        if motioncorr_24params:
            num_confounds=3+24
        else:
            num_confounds=3+6
        WM_signal=extract_mask_trace(self.inputs.bold, self.inputs.WM_mask)
        confounds=np.zeros([np.size(WM_signal,0), num_confounds])
        confounds[:,0]=WM_signal
        confounds[:,1]=extract_mask_trace(self.inputs.bold, self.inputs.CSF_mask)
        confounds[:,2]=extract_mask_trace(self.inputs.bold, self.inputs.brain_mask)
        if motioncorr_24params:
            confounds[:,3:27]=motion_24_params(self.inputs.movpar_file)
            csv_columns=['WM_signal', 'CSF_signal', 'global_signal', 'mov1-1', 'mov2-1', 'mov3-1', 'rot1-1', 'rot2-1', 'rot3-1', 'mov1^2', 'mov2^2', 'mov3^2', 'rot1^2', 'rot2^2', 'rot3^2', 'mov1-1^2', 'mov2-1^2', 'mov3-1^2', 'rot1-1^2', 'rot2-1^2', 'rot3-1^2']

        else:
            confounds[:,3:9]=extract_rigid_movpar(self.inputs.movpar_file)
            csv_columns=['WM_signal', 'CSF_signal', 'global_signal', 'mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3']

        confounds_csv=write_confound_csv(confounds, csv_columns)
        if self.inputs.apply_GSR:
            cleaned=clean_bold(self.inputs.bold, confounds, self.inputs.TR)
        else:
            cleaned=clean_bold(self.inputs.bold, confounds[:, np.r_[0:3,4:]], self.inputs.TR)

        setattr(self, 'cleaned_bold', cleaned)
        setattr(self, 'confounds_csv', confounds_csv)
        return runtime

    def _list_outputs(self):
        return {'cleaned_bold': getattr(self, 'cleaned_bold'),
                'confounds_csv': getattr(self, 'confounds_csv')}

def write_confound_csv(confound_array, column_names):
    import pandas as pd
    import os
    df = pd.DataFrame(confound_array)
    df.columns=column_names
    csv_path=os.path.abspath("confounds.csv")
    df.to_csv(csv_path)
    return csv_path

def clean_bold(bold, confounds_array, TR):
    '''clean with nilearn'''
    import nilearn.image
    import os
    regressed_bold = nilearn.image.clean_img(bold, detrend=True, standardize=True, high_pass=0.01, confounds=confounds_array, t_r=TR)
    cleaned = nilearn.image.smooth_img(regressed_bold, 0.3)
    cleaned_path=os.path.abspath('cleaned.nii.gz')
    cleaned.to_filename(cleaned_path)
    return cleaned_path

def motion_24_params(movpar_csv):
    import numpy as np
    rigid_params=extract_rigid_movpar(movpar_csv)
    movpar=np.zeros([np.size(rigid_params,0), 24])
    movpar[:,:6]=rigid_params
    for i in range(6):
        #add the timepoint 1 TR before
        movpar[0,6+i]=movpar[0,i]
        movpar[1:,6+i]=movpar[:-1,i]
        #add the squared coefficients
        movpar[:,12+i]=movpar[:,i]**2
        movpar[:,18+i]=movpar[:,6+i]**2
    return movpar

def extract_rigid_movpar(movpar_csv):
    import numpy as np
    import csv
    temp = []
    with open(movpar_csv) as csvfile:
        motcorr = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in motcorr:
            temp.append(row)
    movpar=np.zeros([(len(temp)-1), 6])
    j=0
    for row in temp[1:]:
        for i in range(2,len(row)):
            movpar[j,i-2]=float(row[i])
        j=j+1
    return movpar


def extract_mask_trace(bold, mask):
    import numpy as np
    import nilearn.masking
    mask_signal=nilearn.masking.apply_mask(bold, mask)
    mean_trace=np.mean(mask_signal, 1)
    return mean_trace



def extract_labels(atlas):
    import nilearn.regions
    nilearn.regions.connected_label_regions(atlas)


class MaskEPIInputSpec(BaseInterfaceInputSpec):
    mask = File(exists=True, mandatory=True, desc="Mask to transfer to EPI space")
    ref_EPI = File(exists=True, mandatory=True, desc="Ref 3D EPI")
    EPI_to_anat_trans = File(exists=True, desc="Transforms for registration of EPI to anat")
    use_transforms = traits.Bool(mandatory=True, desc="determine whether transform is used")

class MaskEPIOutputSpec(TraitedSpec):
    EPI_mask = traits.File(desc="The generated EPI mask")

class MaskEPI(BaseInterface):

    input_spec = MaskEPIInputSpec
    output_spec = MaskEPIOutputSpec

    def _run_interface(self, runtime):
        import os
        import nibabel as nb
        from nipype.interfaces.base import CommandLine
        resampled_mask_path=os.path.abspath('resampled_mask.nii.gz')
        new_mask_path=os.path.abspath('EPI_mask.nii.gz')
        if self.inputs.use_transforms:
            to_EPI = CommandLine('antsApplyTransforms', args='-i ' + self.inputs.mask + ' -r ' + self.inputs.ref_EPI + ' -t ' + self.inputs.EPI_to_anat_trans + ' -o ' + resampled_mask_path + ' -n GenericLabel')
            to_EPI.run()
        else:
            to_EPI = CommandLine('antsApplyTransforms', args='-i ' + self.inputs.mask + ' -r ' + self.inputs.ref_EPI + ' -o ' + resampled_mask_path + ' -n GenericLabel')
            to_EPI.run()

        nb.Nifti1Image(nb.load(resampled_mask_path).dataobj, nb.load(self.inputs.ref_EPI).affine,
                       nb.load(self.inputs.ref_EPI).header).to_filename(new_mask_path)

        setattr(self, 'EPI_mask', new_mask_path)
        return runtime

    def _list_outputs(self):
        return {'EPI_mask': getattr(self, 'EPI_mask')}
