import sys
#from roseasy import movers
from pyrosetta import *

test = sys.argv[1]

if test=='relax':
    from roseasy.movers import relax as r
    init()
    pdbpath = 'test_inputs/4lhx_chainA.pdb'
    pose = pose_from_file(pdbpath)
    relax = r.Relax()
    relax.rounds = 1
    relax.pose = pose
    relax.apply()

elif test=='fkic':
    from roseasy.movers.loopmodeler import LoopModeler
    init()
    pdbpath = 'test_inputs/4cmp.pdb'
    pose = pose_from_file(pdbpath)
    lm = LoopModeler()
    lm.config = 'fkic'
    lm.frag_files = \
            ['test_inputs/fragments/4cmpA/4cmpA.200.9mers.rewrite.gz',
                    'test_inputs/fragments/4cmpA/4cmpA.200.3mers.rewrite.gz']
    lm.frag_sizes = [3, 9]
    lm.loops_from_file('test_inputs/cas_loops')
    lm.fa_temp_cycles = 1
    lm.cen_temp_cycles = 1
    lm.pose = pose
    print(lm.pose.pdb_info().pose2pdb(710))
    lm.apply()
