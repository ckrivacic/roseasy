from roseasy.movers import add_residues
from roseasy.movers.loopmodeler import LoopModeler
from pyrosetta import init
from pyrosetta import pose_from_file
import os

inputs = '/Users/codykrivacic/cas/roseasy/roseasy/tests/test_inputs/'
init()
pose = pose_from_file(os.path.join(inputs,'4lhx_chainA.pdb'))

add_residues.add_aas(pose, 30, 'VVV', pdbnum=True)
pose.dump_pdb(os.path.join(inputs, '4lhx_added_resis.pdb'))

pose = pose_from_file(os.path.join(inputs, '4lhx_added_resis.pdb'))

lm = LoopModeler()
lm.config='ngk'
lm.mark_as_test_run()
lm.loops_from_file(os.path.join(inputs,'4lhx.loop'))
lm.pose = pose
print('Applying...')
lm.apply()
print("APPLIED")
pose.dump_pdb(os.path.join(inputs, '4lhx_chainA_modeled.pdb'))
