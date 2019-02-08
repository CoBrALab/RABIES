from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


def init_anat_mask_prep_wf(subject_csv,name='anat_prep_mask_wf'):
    '''
    This workflow will take the output masks and labels from pydpyper for each
    subject, the transform of each subject,
    '''

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['csv_file']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_anat']), name='outputnode')

    return workflow
'''
    workflow.connect([
        (inputnode, pydpyper, [("csv_file", "csv_file")]),
        (pydpyper, outputnode, [("out_anat", "out_anat")]),
    ])
'''

'''
path labels:pydpyper/mbm_atlasReg_processed/${var:0:8}preproc_anat/voted.mnc
path mask: pydpyper/mbm_atlasReg_atlases/DSURQE_40micron_mask/resampled/${var:0:8}preproc_anat_I_lsq6_lsq12_and_nlin-resampled_mask.mnc
transform: mbm_atlasReg_processed/MFC_007_preproc_anat/transforms/MFC_007_preproc_anat__concat_lsq6_I_lsq6_lsq12_and_nlin.xfm
'''
