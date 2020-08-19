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
    split_name = traits.Str(mandatory=True, desc="String with info about the main_split.")
    name_source = traits.Str(mandatory=True, desc="Input file template for naming outputs.")


class PlotOverlapOutputSpec(TraitedSpec):
    out_png = File(exists=True, desc="Output png.")


class PlotOverlap(BaseInterface):

    input_spec = PlotOverlapInputSpec
    output_spec = PlotOverlapOutputSpec

    def _run_interface(self, runtime):
        import os
        import pathlib
        folder_template = pathlib.Path(self.inputs.split_name).name.rsplit(".nii")[0]
        filename_template = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")[0]

        import rabies
        dir_path = os.path.dirname(os.path.realpath(rabies.__file__))
        script_path = dir_path+'/shell_scripts/plot_overlap.sh'
        os.makedirs(self.inputs.out_dir+'/'+folder_template, exist_ok=True)
        out_name = self.inputs.out_dir+'/'+folder_template+'/' + \
            filename_template+'_'+self.inputs.reg_name+'.png'

        from rabies.preprocess_pkg.utils import run_command
        command = 'bash %s %s %s %s' % (
            script_path, self.inputs.moving, self.inputs.fixed, out_name)
        rc = run_command(command)

        setattr(self, 'out_png', out_name)
        return runtime

    def _list_outputs(self):
        return {'out_png': getattr(self, 'out_png')}


class PlotMotionTraceInputSpec(BaseInterfaceInputSpec):
    confounds_csv = File(exists=True, mandatory=True,
                         desc="Confound csv file with motion timecourses.")
    out_dir = traits.Str(mandatory=True, desc="Directory for QC outputs.")
    split_name = traits.Str(mandatory=True, desc="String with info about the main_split.")
    name_source = traits.Str(mandatory=True, desc="Input file template for naming outputs.")


class PlotMotionTraceOutputSpec(TraitedSpec):
    out_png = File(exists=True, desc="Output png.")


class PlotMotionTrace(BaseInterface):

    input_spec = PlotMotionTraceInputSpec
    output_spec = PlotMotionTraceOutputSpec

    def _run_interface(self, runtime):
        import os
        import pathlib
        folder_template = pathlib.Path(self.inputs.split_name).name.rsplit(".nii")[0]
        filename_template = pathlib.Path(self.inputs.name_source).name.rsplit(".nii")[0]

        import pandas as pd
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
        os.makedirs(self.inputs.out_dir+'/'+folder_template, exist_ok=True)
        prefix = self.inputs.out_dir+'/'+folder_template+'/' + \
            filename_template
        command = 'bash %s %s %s' % (script_path, par_file, prefix)
        from rabies.preprocess_pkg.utils import run_command
        rc = run_command(command)

        setattr(self, 'out_png', '%s_motion_traces.png' % (prefix))
        return runtime

    def _list_outputs(self):
        return {'out_png': getattr(self, 'out_png')}
