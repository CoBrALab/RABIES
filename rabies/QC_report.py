from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, BaseInterface
)


class PlotOverlapInputSpec(BaseInterfaceInputSpec):
    moving = File(exists=True, mandatory=True,
                  desc="Moving image from registration.")
    fixed = File(exists=True, mandatory=True,
                 desc="Fixed image from registration.")
    reg_name = traits.Str(
        mandatory=True, desc="Name of the registration which is displayed.")
    out_dir = traits.Str(mandatory=True, desc="Directory for QC outputs.")


class PlotOverlapOutputSpec(TraitedSpec):
    out_png = File(exists=True, desc="Output png.")


class PlotOverlap(BaseInterface):

    input_spec = PlotOverlapInputSpec
    output_spec = PlotOverlapOutputSpec

    def _run_interface(self, runtime):
        import os
        # if not os.environ["template_anat"]==self.inputs.fixed:
        if not 'resampled_template.nii.gz' == os.path.basename(self.inputs.fixed):
            subject_id = 'sub-' + \
                (os.path.basename(self.inputs.moving).split(
                    '_ses-')[0]).split('sub-')[1]
            session = os.path.basename(self.inputs.moving).split('_ses-')[1][0]
            if 'run' in os.path.basename(self.inputs.moving):
                run = os.path.basename(self.inputs.moving).split('_run-')[1][0]
                filename_template = '%s_ses-%s_run-%s' % (
                    subject_id, session, run)
            else:
                filename_template = '%s_ses-%s' % (subject_id, session)

        import rabies
        dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
        script_path = dir_path+'/shell_scripts/plot_overlap.sh'
        # if not os.environ["template_anat"]==self.inputs.fixed:
        if not 'resampled_template.nii.gz' == os.path.basename(self.inputs.fixed):
            os.makedirs(self.inputs.out_dir+'/'+subject_id, exist_ok=True)
            out_name = self.inputs.out_dir+'/'+subject_id+'/' + \
                filename_template+'_'+self.inputs.reg_name+'.png'
        else:
            os.makedirs(self.inputs.out_dir, exist_ok=True)
            out_name = self.inputs.out_dir+'/'+self.inputs.reg_name+'.png'
        command = 'bash %s %s %s %s' % (
            script_path, self.inputs.moving, self.inputs.fixed, out_name)
        from rabies.preprocess_pkg.utils import run_command
        rc = run_command(command)

        setattr(self, 'out_png', out_name)
        return runtime

    def _list_outputs(self):
        return {'out_png': getattr(self, 'out_png')}


class PlotMotionTraceInputSpec(BaseInterfaceInputSpec):
    confounds_csv = File(exists=True, mandatory=True,
                         desc="Confound csv file with motion timecourses.")
    out_dir = traits.Str(mandatory=True, desc="Directory for QC outputs.")


class PlotMotionTraceOutputSpec(TraitedSpec):
    out_png = File(exists=True, desc="Output png.")


class PlotMotionTrace(BaseInterface):

    input_spec = PlotMotionTraceInputSpec
    output_spec = PlotMotionTraceOutputSpec

    def _run_interface(self, runtime):
        import os
        import pandas as pd
        subject_id = os.path.basename(
            self.inputs.confounds_csv).split('_ses-')[0]
        session = os.path.basename(
            self.inputs.confounds_csv).split('_ses-')[1][0]
        run = os.path.basename(self.inputs.confounds_csv).split('_run-')[1][0]
        filename_template = '%s_ses-%s_run-%s' % (subject_id, session, run)

        def csv2par(in_confounds):
            df = pd.read_csv(in_confounds)
            new_df = pd.DataFrame(
                columns=['mov1', 'mov2', 'mov3', 'rot1', 'rot2', 'rot3'])
            new_df['mov1'] = df['mov1']
            new_df['mov2'] = df['mov2']
            new_df['mov3'] = df['mov3']
            new_df['rot1'] = df['rot1']
            new_df['rot2'] = df['rot2']
            new_df['rot3'] = df['rot3']
            out_confounds = os.path.abspath(
                (os.path.basename(in_confounds).split('.')[0])+('.par'))
            new_df.to_csv(out_confounds, sep='\t', index=False, header=False)
            return out_confounds
        par_file = csv2par(self.inputs.confounds_csv)

        import rabies
        dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
        script_path = dir_path+'/shell_scripts/plot_motion_traces.sh'
        os.makedirs(self.inputs.out_dir+'/'+subject_id, exist_ok=True)
        prefix = self.inputs.out_dir+'/'+subject_id+'/'+filename_template
        command = 'bash %s %s %s' % (script_path, par_file, prefix)
        from rabies.preprocess_pkg.utils import run_command
        rc = run_command(command)

        setattr(self, 'out_png', '%s_motion_traces.png' % (prefix))
        return runtime

    def _list_outputs(self):
        return {'out_png': getattr(self, 'out_png')}
