'''
Find RMSD of all files in a folder and compare to a native or input
structure.

Usage:
    rmsd.py <folder> <input>

To do: Add option for recursive lookup.
'''

import docopt
from pyrosetta.rosetta.core.pose import setPoseExtraScore
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta import pose_from_file
from pyrosetta import init
from roseasy import pipeline
import os, gzip

def add_lines_to_gzip(fname, lines):
    final_lines = []
    with gzip.open(fname, 'rt') as f:
        for line in f:
            final_lines.append(line)

    for line in lines:
        final_lines.append(line + '\n')

    with gzip.open(fname, 'wt') as f:
        for line in final_lines:
            f.write(line)

def main():
    args = docopt.docopt(__doc__)
    folder = args['<folder>']
    init()
    target = pose_from_file(args['<input>'])
    workspace = pipeline.workspace_from_dir(folder)
    loop = workspace.largest_loop

    for root, dirs, files in os.walk(workspace.output_dir):
        for name in files:
            if name.endswith('.pdb.gz'):
                pdbpath = os.path.join(root, name)
                mobile = pose_from_file(pdbpath)
                rmsd = CA_rmsd(mobile, target, int(loop.start), int(loop.end))

                metric_name = 'EXTRA_METRIC_Loop_CA_RMSD [[-]]'

                add_lines_to_gzip(pdbpath, [metric_name + ' ' + str(rmsd)])
                #rmsd = all_atom_rmsd(mobile, target)
                #print(rmsd)


if __name__=='__main__':
    main()
