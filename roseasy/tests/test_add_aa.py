from roseasy.movers.add_residues import *

init()
pose = pose_from_file('test_inputs/4cmp.pdb')
add_aas(pose, 1085, 'PPPPP', pdbnum=True)

pose.dump_pdb('test_inputs/prolines.pdb')
