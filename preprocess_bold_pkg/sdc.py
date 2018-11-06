import os
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import Function
from .utils import init_bold_reference_wf, Merge

def init_sdc_wf(tool, name='sdc_wf'):

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['ref_bold', 'reversed_bold_file', 'name_source']),
                        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['fieldwarp', 'affine_trans', 'nlin_trans', 'avg_corrected_image', 'ref_reversed_bold_file']),
        name='outputnode')

    bold_reference_wf = init_bold_reference_wf(enhance_contrast=True)
    workflow.connect([
        (inputnode, bold_reference_wf, [('reversed_bold_file', 'inputnode.bold_file')]),
        (bold_reference_wf, outputnode, [('outputnode.enhanced_ref_image', 'ref_reversed_bold_file')]),
    ])

    if tool=='FSL':
        mk_list = pe.Node(Function(input_names=["first_img", "second_img"],
                       output_names=["list"],
                       function=makeList), name='mk_list')

        from .utils import Merge
        merge = pe.Node(Merge(), name='merge')

        from nipype.interfaces.fsl import TOPUP
        sdc_topup = pe.Node(TOPUP(config=os.path.abspath("/topup/empty.cnf"), encoding_file = os.path.abspath("/topup/topup_encoding.txt"),
            output_type = "NIFTI_GZ", warp_res = 1, fwhm=0.3, max_iter=10, ssqlambda=1, regmod='bending_energy', estmov=1, minmet=0,
            splineorder=3, numprec='double', interp='spline', scale=1
            #topup.inputs.subsamp = 2, topup.inputs.reg_lambda = 0.005
            ), name='SDC_TOPUP')

        workflow.connect([
            (inputnode, mk_list, [('ref_bold', 'first_img')]),
            (bold_reference_wf, mk_list, [('outputnode.ref_image', 'second_img')]),
            (mk_list, merge, [('list', 'in_files')]),
            (inputnode, merge, [('name_source', 'header_source')]),
            (merge, sdc_topup, [('out_file', 'in_file')]),
            (sdc_topup, outputnode, [
                ('out_field', 'fieldwarp'),
                ('out_corrected', 'avg_corrected_image')]),
        ])

    elif tool=='ANTs':
        from .utils import antsGenerateTemplate
        sdc_antsGenTemplate = pe.Node(antsGenerateTemplate(), name='sdc_antsGenTemplate')
        workflow.connect([
            (inputnode, sdc_antsGenTemplate, [('ref_bold', 'EPI')]),
            (bold_reference_wf, sdc_antsGenTemplate, [('outputnode.ref_image', 'reversed_EPI')]),
            (sdc_antsGenTemplate, outputnode, [
                ('affine_trans', 'affine_trans'),
                ('nlin_trans', 'nlin_trans')]),
        ])


    return workflow


def writeCSV(first_img, second_img):
    #create a CSV file as input for antsGenerateTemplate
    import csv
    import os
    with open('EPI.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerow([first_img])
        csvwriter.writerow([second_img])
    inputs = os.path.abspath('EPI.csv')

    return inputs

def makeList(first_img, second_img):
    #merge two files into a list
    file_list = [first_img, second_img]
    return file_list
