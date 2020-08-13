from roseasy.standard_params import filters
import os
from pyrosetta import pose_from_file
from pyrosetta import init
from roseasy.pipeline import workspace_from_dir


rosetta_dir = '~/software/rosetta'
ws = workspace_from_dir(os.path.expanduser('~/cas/test/'))
dalphaball_path = os.path.join(rosetta_dir, 'source',
        'external', 'DAlphaBall', 'DAlphaBall.gcc')
init('-holes:dalphaball {}'.format(dalphaball_path))
pose =  pose_from_file('/Users/codykrivacic/cas/roseasy/roseasy/tests/test_inputs/4lhx_chainA.pdb')

filters = filters.get_filters(ws)
for f in filters:
    f.apply(pose)
