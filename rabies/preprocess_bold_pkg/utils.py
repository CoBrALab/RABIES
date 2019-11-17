import os
import nibabel as nb
import numpy as np
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, InputMultiPath, BaseInterface, SimpleInterface
)
from nipype.interfaces.base import CommandLine, CommandLineInputSpec

def prep_bids_iter(layout):
    subject_list=layout.get_subject()
    #create a dictionary with list of bold session and run numbers for each subject
    session_iter={}
    run_iter={}
    for sub in subject_list:
        sub_func=layout.get(datatype='func', extension=['nii', 'nii.gz'])
        session=0
        run=0
        for func_bids in sub_func:
            if int(func_bids.get_entities()['session'])>session:
                session=int(func_bids.get_entities()['session'])
            if int(func_bids.get_entities()['run'])>run:
                run=int(func_bids.get_entities()['run'])
        session_iter[sub] = list(range(1,int(session)+1))
        run_iter[sub] = list(range(1,int(run)+1))
    return subject_list, session_iter, run_iter

class BIDSDataGraberInputSpec(BaseInterfaceInputSpec):
    bids_dir = traits.Str(exists=True, mandatory=True, desc="BIDS data directory")
    datatype = traits.Str(exists=True, mandatory=True, desc="datatype of the target file")
    subject_id = traits.Str(exists=True, mandatory=True, desc="Subject ID")
    session = traits.Int(exists=True, mandatory=True, desc="Session number")
    run = traits.Int(exists=True, default=None, desc="Run number")

class BIDSDataGraberOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="Selected file based on the provided parameters.")

class BIDSDataGraber(BaseInterface):
    """
    """

    input_spec = BIDSDataGraberInputSpec
    output_spec = BIDSDataGraberOutputSpec

    def _run_interface(self, runtime):
        import os
        from bids.layout import BIDSLayout
        layout = BIDSLayout(self.inputs.bids_dir)
        try:
            if self.inputs.datatype=='func':
                bids_file=layout.get(subject=self.inputs.subject_id, session=self.inputs.session, run=self.inputs.run, extension=['nii', 'nii.gz'], datatype=self.inputs.datatype)
                func=layout.get(subject=self.inputs.subject_id, session=self.inputs.session, run=self.inputs.run, extension=['nii', 'nii.gz'], datatype=self.inputs.datatype, return_type='filename')
                file=func[0]
            elif self.inputs.datatype=='anat':
                bids_file=layout.get(subject=self.inputs.subject_id, session=self.inputs.session, extension=['nii', 'nii.gz'], datatype=self.inputs.datatype)
                anat=layout.get(subject=self.inputs.subject_id, session=self.inputs.session, extension=['nii', 'nii.gz'], datatype=self.inputs.datatype, return_type='filename')
                file=anat[0]
            else:
                raise ValueError('Wrong datatype %s' % (self.inputs.datatype))
            if len(bids_file)>1:
                raise ValueError('Provided BIDS spec lead to duplicates: %s' % (str(self.inputs.datatype+'_'+self.inputs.subject_id+'_'+self.inputs.session+'_'+self.inputs.run)))
        except:
            raise ValueError('Error with BIDS spec: %s' % (str(self.inputs.datatype+'_'+self.inputs.subject_id+'_'+self.inputs.session+'_'+self.inputs.run)))

        nii_format=bids_file[0].get_entities()['extension']
        #RABIES only work with compressed .nii for now
        if nii_format=='nii':
            os.system('gzip %s' % (file,))
            file=file+'.gz'

        setattr(self, 'out_file', file)

        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}


def init_bold_reference_wf(name='gen_bold_ref'):
    """
    This workflow generates reference BOLD images for a series

    **Parameters**

        bold_file : str
            BOLD series NIfTI file
        name : str
            Name of workflow (default: ``bold_reference_wf``)

    **Inputs**

        bold_file
            BOLD series NIfTI file

    **Outputs**

        bold_file
            Validated BOLD series NIfTI file
        ref_image
            Reference image to which BOLD series is motion corrected
        enhanced_ref_image
            Reference image with enhanced contrast for EPI to anat coregistration
        skip_vols
            Number of non-steady-state volumes detected at beginning of ``bold_file``
        validation_report
            HTML reportlet indicating whether ``bold_file`` had a valid affine


    **Subworkflows**

        * :py:func:`~fmriprep.workflows.bold.util.init_enhance_and_skullstrip_wf`

    """
    from nipype.pipeline import engine as pe
    from nipype.interfaces import utility as niu
    from .interfaces import ValidateImage


    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['bold_file']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['bold_file', 'skip_vols', 'ref_image', 'validation_report', 'enhanced_ref_image']),
        name='outputnode')


    '''
    Check the correctness of x-form headers (matrix and code)

        This interface implements the `following logic
        <https://github.com/poldracklab/fmriprep/issues/873#issuecomment-349394544>
    '''
    validate = pe.Node(ValidateImage(), name='validate')
    validate.plugin_args = {'qsub_args': '-pe smp %s' % (str(2*int(os.environ["min_proc"]))), 'overwrite': True}

    gen_ref = pe.Node(EstimateReferenceImage(), name='gen_ref')
    gen_ref.plugin_args = {'qsub_args': '-pe smp %s' % (str(2*int(os.environ["min_proc"]))), 'overwrite': True}

    workflow.connect([
        (inputnode, validate, [('bold_file', 'in_file')]),
        (validate, gen_ref, [('out_file', 'in_file')]),
        (validate, outputnode, [('out_file', 'bold_file'),
                                ('out_report', 'validation_report')]),
        (gen_ref, outputnode, [('ref_image', 'ref_image'),
                               ('n_volumes_to_discard', 'skip_vols')]),
    ])

    return workflow


class EstimateReferenceImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="4D EPI file")

class EstimateReferenceImageOutputSpec(TraitedSpec):
    ref_image = File(exists=True, desc="3D reference image")
    n_volumes_to_discard = traits.Int(desc="Number of detected non-steady "
                                           "state volumes in the beginning of "
                                           "the input file")

class EstimateReferenceImage(BaseInterface):
    """
    Given an 4D EPI file estimate an optimal reference image that could be later
    used for motion estimation and coregistration purposes. If detected uses
    anat saturated volumes (non-steady state). Otherwise, a median of
    of a subset of motion corrected volumes is used. In the later case, a first
    median is extracted from the raw data and used as reference for motion correction,
    then a new median image is extracted from the corrected series, and the process
    is repeated one more time to generate a final image reference image.
    """

    input_spec = EstimateReferenceImageInputSpec
    output_spec = EstimateReferenceImageOutputSpec

    def _run_interface(self, runtime):

        import os
        import nibabel as nb
        import numpy as np

        in_nii = nb.load(self.inputs.in_file)
        data_slice = in_nii.dataobj[:, :, :, :50]

        # Slicing may induce inconsistencies with shape-dependent values in extensions.
        # For now, remove all. If this turns out to be a mistake, we can select extensions
        # that don't break pipeline stages.
        in_nii.header.extensions.clear()

        n_volumes_to_discard = _get_vols_to_discard(in_nii)

        subject_id=os.path.basename(self.inputs.in_file).split('_ses-')[0]
        session=os.path.basename(self.inputs.in_file).split('_ses-')[1][0]
        run=os.path.basename(self.inputs.in_file).split('_run-')[1][0]
        filename_template = '%s_ses-%s_run-%s' % (subject_id, session, run)

        out_ref_fname = os.path.abspath('%s_bold_ref.nii.gz' % (filename_template))

        if n_volumes_to_discard == 0:
            print("Detected no dummy scans. Generating the ref EPI based on multiple volumes.")
            #if no dummy scans, will generate a median from a subset of max 40
            #slices of the time series
            if in_nii.shape[-1] > 40:
                slice_fname = os.path.abspath("slice.nii.gz")
                nb.Nifti1Image(data_slice[:, :, :, 20:40], in_nii.affine,
                               in_nii.header).to_filename(slice_fname)
                median_fname = os.path.abspath("median.nii.gz")
                nb.Nifti1Image(np.median(data_slice[:, :, :, 20:40], axis=3), in_nii.affine,
                               in_nii.header).to_filename(median_fname)
            else:
                slice_fname = self.inputs.in_file
                median_fname = os.path.abspath("median.nii.gz")
                nb.Nifti1Image(np.median(data_slice, axis=3), in_nii.affine,
                               in_nii.header).to_filename(median_fname)

            print("First iteration to generate reference image.")
            res = antsMotionCorr(in_file=slice_fname, ref_file=median_fname, second=False).run()
            median = np.median(nb.load(res.outputs.mc_corrected_bold).get_data(), axis=3)
            tmp_median_fname = os.path.abspath("tmp_median.nii.gz")
            nb.Nifti1Image(median, in_nii.affine,
                           in_nii.header).to_filename(tmp_median_fname)

            print("Second iteration to generate reference image.")
            res = antsMotionCorr(in_file=slice_fname, ref_file=tmp_median_fname, second=True).run()
            median_image_data = np.median(nb.load(res.outputs.mc_corrected_bold).get_data(), axis=3)
        else:
            print("Detected "+str(n_volumes_to_discard)+" dummy scans. Taking the median of these volumes as reference EPI.")
            median_image_data = np.median(
                data_slice[:, :, :, :n_volumes_to_discard], axis=3)

        #median_image_data is a 3D array of the median image, so creates a new nii image
        #saves it
        nb.Nifti1Image(median_image_data, in_nii.affine,
                       in_nii.header).to_filename(out_ref_fname)


        setattr(self, 'ref_image', out_ref_fname)
        setattr(self, 'n_volumes_to_discard', n_volumes_to_discard)

        return runtime

    def _list_outputs(self):
        return {'ref_image': getattr(self, 'ref_image'),
                'n_volumes_to_discard': getattr(self, 'n_volumes_to_discard')}


'''
Takes a nifti file, extracts the mean signal of the first 50 volumes and computes which are outliers.
is_outlier function: computes Modified Z-Scores (https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm) to determine which volumes are outliers.
'''
def _get_vols_to_discard(img):
    from nipype.algorithms.confounds import is_outlier
    data_slice = img.dataobj[:, :, :, :50]
    global_signal = data_slice.mean(axis=0).mean(axis=0).mean(axis=0)
    return is_outlier(global_signal)


class antsMotionCorrInputSpec(CommandLineInputSpec):
    in_file = File(exists=True, mandatory=True, argstr='%s, 1 , 20, Regular, 0.2] -t Rigid[0.25] -i 50x20 -u 1 -e 1 -s 1x0 -f 2x1 -n 10 -l 1 -w 1', position=2, desc='input BOLD time series')
    ref_file = File(exists=True, mandatory=True, argstr='-d 3 -o [ ants_mc_tmp/motcorr, ants_mc_tmp/motcorr.nii.gz, ants_mc_tmp/motcorr_avg.nii.gz] -m mi[ %s, ', position=1, desc='ref file to realignment time series')
    second = traits.Bool(desc="specify if it is the second iteration")

class antsMotionCorrOutputSpec(TraitedSpec):
    mc_corrected_bold = File(exists=True, desc="motion corrected time series")
    motcorr_warp = File(exists=True, desc="motion estimation of the time series")
    avg_image = File(exists=True, desc="average image of the motion corrected time series")
    csv_params = File(exists=True, desc="csv files with the 6-parameters rigid body transformations")

class antsMotionCorr(CommandLine):
    _cmd = 'antsMotionCorr'
    input_spec = antsMotionCorrInputSpec
    output_spec = antsMotionCorrOutputSpec

    def _run_interface(self, runtime):

        #change the name of the first iteration directory to prevent overlap of files with second iteration
        if self.inputs.second:
            mk_tmp = CommandLine('mv', args='ants_mc_tmp first_ants_mc_tmp')
            mk_tmp.run()

        #make a tmp directory to store the files
        import os
        os.makedirs('ants_mc_tmp', exist_ok=True)

        # Run the command line as a natural CommandLine interface
        runtime = super(antsMotionCorr, self)._run_interface(runtime)

        setattr(self, 'csv_params', 'ants_mc_tmp/motcorrMOCOparams.csv')
        setattr(self, 'motcorr_warp', 'ants_mc_tmp/motcorrWarp.nii.gz')
        setattr(self, 'mc_corrected_bold', 'ants_mc_tmp/motcorr.nii.gz')
        setattr(self, 'avg_image', 'ants_mc_tmp/motcorr_avg.nii.gz')

        return runtime

    def _list_outputs(self):
        return {'mc_corrected_bold': getattr(self, 'mc_corrected_bold'),
                'motcorr_warp': getattr(self, 'motcorr_warp'),
                'csv_params': getattr(self, 'csv_params'),
                'avg_image': getattr(self, 'avg_image')}


class antsGenerateTemplateInputSpec(CommandLineInputSpec):
    #optional further info -f 4x2x1 -s 2x1x0vox -q 30x20x4
    EPI = File(exists=True, mandatory=True, argstr='-d 3 -r 1 -b 1 -c 0 -k 1 -t SyN -m CC -o sdc_ %s', position=0, desc='reference image from the original EPI')
    reversed_EPI = File(exists=True, mandatory=True, argstr=' %s', position=1, desc='generated reference image from the reversed EPI')

class antsGenerateTemplateOutputSpec(TraitedSpec):
    affine_trans = File(exists=True, desc="affine transformations to the template")
    nlin_trans = File(exists=True, desc="non-linear transforms to the template")
    template = File(exists=True, desc="generated template file")

class antsGenerateTemplate(CommandLine):
    _cmd = 'antsMultivariateTemplateConstruction2.sh'
    input_spec = antsGenerateTemplateInputSpec
    output_spec = antsGenerateTemplateOutputSpec

    def _run_interface(self, runtime):

        # Run the command line as a natural CommandLine interface
        runtime = super(antsGenerateTemplate, self)._run_interface(runtime)

        setattr(self, 'affine_trans', os.path.abspath("sdc_ref_EPI00GenericAffine.mat"))
        setattr(self, 'nlin_trans', os.path.abspath("sdc_ref_EPI01Warp.nii.gz"))
        setattr(self, 'template', os.path.abspath("sdc_template0.nii.gz"))

        return runtime

    def _list_outputs(self):
        return {'affine_trans': getattr(self, 'affine_trans'),
                'nlin_trans': getattr(self, 'nlin_trans'),
                'template': getattr(self, 'template')}


class slice_applyTransformsInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="Input 4D EPI")
    ref_file = File(exists=True, mandatory=True, desc="The reference 3D space to which the EPI will be warped.")
    transforms = traits.List(desc="List of transforms to apply to every slice")
    inverses = traits.List(desc="Define whether some transforms must be inverse, with a boolean list where true defines inverse e.g.[0,1,0]")
    apply_motcorr = traits.Bool(default=True, desc="Whether to apply motion realignment.")
    motcorr_params = File(exists=True, desc="xforms from head motion estimation .csv file")
    isotropic_resampling = traits.Bool(desc="If true, the EPI will be resampled to isotropic resolution based on the lowest dimension.")
    upsampling = traits.Float(default=1.0, desc="Option to upsample the voxel resolution upon resampling to minimize data loss. All dimensions will be multiplied by the specified proportion. e.g. 2.0 doubles the resolution.")
    data_type = traits.Str(default='float64', desc="Specify resampling data format to control for file size. Can specify a numpy data type from https://docs.scipy.org/doc/numpy/user/basics.types.html.")

class slice_applyTransformsOutputSpec(TraitedSpec):
    out_files = traits.List(desc="warped images after the application of the transforms")

class slice_applyTransforms(BaseInterface):
    """
    This interface will apply head motion correction as well as susceptibility distortion correction
    if specified to the input EPI volumes using antsApplyTransforms.
    """

    input_spec = slice_applyTransformsInputSpec
    output_spec = slice_applyTransformsOutputSpec

    def _run_interface(self, runtime):
        #resampling the reference image to the dimension of the EPI
        from nibabel import processing
        import numpy as np
        import nibabel as nb
        import os
        img=nb.load(self.inputs.in_file)
        shape=img.header.get_zooms()[:3]
        upsampling_multiplier=1/self.inputs.upsampling
        if self.inputs.isotropic_resampling:
            low_dim=np.asarray(shape).min()
            shape=(low_dim,low_dim,low_dim)
        processing.resample_to_output(nb.load(self.inputs.ref_file), voxel_sizes=(shape[0]*upsampling_multiplier,shape[1]*upsampling_multiplier,shape[2]*upsampling_multiplier), order=4).to_filename('resampled.nii.gz')

        #tranforms is a list of transform files, set in order of call within antsApplyTransforms
        transform_string=""
        for transform,inverse in zip(self.inputs.transforms, self.inputs.inverses):
            if bool(inverse):
                transform_string += "-t [%s,1] " % (transform,)
            else:
                transform_string += "-t %s " % (transform,)

        print("Splitting bold file into lists of single volumes")
        [bold_volumes, num_volumes] = split_volumes(self.inputs.in_file, "bold_")

        if self.inputs.apply_motcorr:
            motcorr_params=self.inputs.motcorr_params
        ref_img=os.path.abspath('resampled.nii.gz')
        warped_volumes = []
        for x in range(0, num_volumes):
            warped_vol_fname = os.path.abspath("deformed_volume" + str(x) + ".nii.gz")
            warped_volumes.append(warped_vol_fname)
            if self.inputs.apply_motcorr:
                os.system('antsMotionCorrStats -m %s -o motcorr_vol%s.mat -t %s' % (motcorr_params, x, x))
                os.system('antsApplyTransforms -i %s %s-t motcorr_vol%s.mat -n BSpline[5] -r %s -o %s' % (bold_volumes[x], transform_string, x, ref_img, warped_vol_fname))
            else:
                os.system('antsApplyTransforms -i %s %s-n BSpline[5] -r %s -o %s' % (bold_volumes[x], transform_string, ref_img, warped_vol_fname))
            #change image to specified data type
            img=nb.load(warped_vol_fname)
            img.set_data_dtype(self.inputs.data_type)
            nb.save(img, warped_vol_fname)

        setattr(self, 'out_files', warped_volumes)
        return runtime

    def _list_outputs(self):
        return {'out_files': getattr(self, 'out_files')}

def split_volumes(in_file, output_prefix):
    '''
    Takes as input a 4D .nii file and splits it into separate time series
    volumes by splitting on the 4th dimension
    '''
    import os
    import numpy as np
    import nibabel as nb
    in_nii = nb.load(in_file)
    num_dimensions = len(in_nii.shape)
    num_volumes = in_nii.shape[3]

    if num_dimensions!=4:
        print("the input file must be of dimensions 4")
        return None

    volumes = []
    for x in range(0, num_volumes):
        data_slice = in_nii.dataobj[:, :, :, x]
        slice_fname = os.path.abspath(output_prefix + "vol" + str(x) + ".nii.gz")
        nb.Nifti1Image(data_slice, in_nii.affine,
                       in_nii.header).to_filename(slice_fname)
        volumes.append(slice_fname)

    return [volumes, num_volumes]



class MergeInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath(File(exists=True), mandatory=True,
                              desc='input list of files to merge, listed in the order to merge')
    header_source = File(exists=True, mandatory=True, desc='a Nifti file from which the header should be copied')
    data_type = traits.Str(default='float64', desc="Specify resampling data format to control for file size. Can specify a numpy data type from https://docs.scipy.org/doc/numpy/user/basics.types.html.")

class MergeOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output merged file')

class Merge(BaseInterface):
    """
    Takes a list of 3D Nifti files and merge them in the order listed.
    """

    input_spec = MergeInputSpec
    output_spec = MergeOutputSpec

    def _run_interface(self, runtime):
        import os
        import nibabel as nb
        import numpy as np

        subject_id=os.path.basename(self.inputs.header_source).split('_ses-')[0]
        session=os.path.basename(self.inputs.header_source).split('_ses-')[1][0]
        run=os.path.basename(self.inputs.header_source).split('_run-')[1][0]
        filename_template = '%s_ses-%s_run-%s' % (subject_id, session, run)

        img = nb.load(self.inputs.in_files[0])
        affine = img.affine
        header = nb.load(self.inputs.header_source).header
        length = len(self.inputs.in_files)
        combined = np.zeros((img.shape[0], img.shape[1], img.shape[2], length))

        i=0
        for file in self.inputs.in_files:
            combined[:,:,:,i] = nb.load(file).dataobj[:,:,:]
            i = i+1
        if (i!=length):
            print("Error occured with Merge.")
            return None
        combined_files = os.path.abspath("%s_combined.nii.gz" % (filename_template))
        combined_image=nb.Nifti1Image(combined, affine,
                       header)
        #change image to specified data type
        combined_image.set_data_dtype(self.inputs.data_type)
        nb.save(combined_image, combined_files)

        setattr(self, 'out_file', combined_files)
        return runtime

    def _list_outputs(self):
        return {'out_file': getattr(self, 'out_file')}



class SkullstripInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="4D EPI file")
    brain_mask = File(exists=True, mandatory=True, desc="Brain mask for brain extraction")

class SkullstripOutputSpec(TraitedSpec):
    skullstrip_brain = File(exists=True, desc="Extracted brain")


class Skullstrip(BaseInterface):

    input_spec = SkullstripInputSpec
    output_spec = SkullstripOutputSpec

    def _run_interface(self, runtime):

        import os
        import nibabel as nb
        file_path=os.path.abspath("skullstrip.nii.gz")
        from nilearn.input_data import NiftiMasker
        nifti_masker = NiftiMasker(
            detrend=False,
            standardize=False,
            mask_img=self.inputs.brain_mask,
            memory='nilearn_cache', memory_level=1)  # cache options
        masked = nifti_masker.fit_transform(self.inputs.in_file)
        nifti_masker.inverse_transform(masked).to_filename(file_path)

        setattr(self, 'skullstrip_brain', file_path)

        return runtime

    def _list_outputs(self):
        return {'skullstrip_brain': getattr(self, 'skullstrip_brain')}
