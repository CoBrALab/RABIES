from bids.layout import BIDSLayout
layout = BIDSLayout("/data/ds000114/")


from nipype.interfaces.io import BIDSDataGrabber
from nipype.pipeline import Node, MapNode, Workflow
from nipype.interfaces.utility import Function


bg_all = Node(BIDSDataGrabber(), name='bids-grabber')
bg_all.inputs.base_dir = '/data/ds000114' #give the path to BIDS dataset directory
bg_all.inputs.output_query = {'bolds': dict(type='bold')}
bg_all.iterables = ('subject', layout.get_subjects()[:2])
wf = Workflow(name="bids_demo")
wf.connect(bg_all, "bolds", analyzeBOLD, "paths")
wf.run()
