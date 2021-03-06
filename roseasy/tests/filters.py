from roseasy.standard_params import filters
import os, sys
from pyrosetta import pose_from_file
from pyrosetta import init
from roseasy.pipeline import workspace_from_dir


#rosetta_dir = '~/software/rosetta'
ws = workspace_from_dir(os.path.expanduser('~/cas/test/'))
print(ws.largest_loop)
print(ws.loops_path)
#sys.exit()
dalphaball_path = os.path.join(ws.rosetta_dir, 'source',
        'external', 'DAlphaBall', 'DAlphaBall.gcc')
init('-holes:dalphaball {} -in:file:s {}'.format(dalphaball_path,
    ws.input_pdb_path))
pose =  pose_from_file(ws.input_pdb_path)

filters = filters.get_filters(ws, score_fragments=True, test_run=True)
#for f in filters:
filters[-1].apply(pose)

pose.dump_pdb('test_out.pdb')
