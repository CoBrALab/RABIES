from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from os.path import join as opj
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.interfaces.utility import Function


def init_bold_iter_wf(data_dir_path, name='init_bold_iter_wf'):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['subject_id', 'session', 'run']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['bold_file']), name='outputnode')

    run_iter = pe.Node(niu.IdentityInterface(fields=['run']),
                      name="run_iter")
    run_num=1
    run_iter.iterables = [('run', list(range(1,run_num+1)))]

    bold_file = opj('{subject_id}', 'ses-{session_num}', 'bold', '{subject_id}_ses-{session_num}_run-{run_num}_bold.nii.gz')

    bold_selectfiles = pe.Node(SelectFiles({'bold': bold_file},
                                   base_directory=data_dir_path),
                       name="bold_selectfiles")


    workflow.connect([
        (inputnode, bold_selectfiles, [
            ("subject_id", "subject_id"),
            ("session", "session_num"),
            ]),
        (run_iter, bold_selectfiles, [
            ("run", "run_num"),
            ]),
        (bold_selectfiles, outputnode, [("bold", "bold_file")]),
    ])

    return workflow
