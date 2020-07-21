import sys
#from roseasy import movers
from pyrosetta import *

test = sys.argv[1]
print('HDIOFSJH:LDSK J:FAUDISLKFC JLKSDA')

if test=='relax':
    from roseasy.movers import relax as r
    init()
    pdbpath = 'test_inputs/4lhx_chainA.pdb'
    pose = pose_from_file(pdbpath)
    relax = r.Relax()
    relax.rounds = 1
    relax.pose = pose
    relax.apply()
