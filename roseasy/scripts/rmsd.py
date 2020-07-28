'''
Find RMSD of all files in a folder and compare to a native or input
structure.

Usage:
    rmsd.py <folder> <input>

To do: Add option for recursive lookup.
'''

import docopt
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta import pose_from_file
from pyrosetta import init
import os

def main():
    args = docopt.docopt(__doc__)
    folder = args['<folder>']
    init()
    target = pose_from_file(args['<input>'])

    for root, dirs, files in os.walk(folder):
        for name in files:
            if name.endswith('.pdb'):
                mobile = pose_from_file(os.path.join(root, name))
                rmsd = CA_rmsd(mobile, target)
                print(rmsd)
                rmsd = all_atom_rmsd(mobile, target)
                print(rmsd)


if __name__=='__main__':
    main()
