from roseasy.scripts import filters
import os
from pyrosetta import pose_from_file
from pyrosetta import init


rosetta_dir = '~/software/rosetta'
dalphaball_path = os.path.join(rosetta_dir, 'source',
        'external', 'DAlphaBall', 'DAlphaBall.gcc')
init('-holes:dalphaball {}'.format(dalphaball_path))
pose =  pose_from_file('/Users/codykrivacic/cas/roseasy/roseasy/tests/test_inputs/4lhx_chainA.pdb')

filters = filters.get_filters()
for f in filters:
    f.apply(pose)
