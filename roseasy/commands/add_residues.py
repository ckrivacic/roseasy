"""
Insert residues into a pdb file and do a quick minimization to attempt
loop closure. You may pass either a single PDB file or a directory of
PDB files. 
Warning: This will move your input PDB files to a backup
outpur directory if a directory is given, or to '<filename>.pdb.bk' if
a single file is given, and overwrite the original filename with the
inserted loop.
For now, this only works on a single folder. In the future it may be
useful to add the capability to pass any number of folders as input.

Usage:
    roseasy add_residues <file_or_folder> <residue_string> <pdb_position> [options]

Options:
    --chain [CH], -c [CH]:
        Specify a chain for the insertion. Defaults to A.  [default: A]
    --out [STR], -o [STR]
        Specify an output file

"""

import os, shutil, glob
import docopt
from pyrosetta import init
from pyrosetta import pose_from_file
from roseasy.movers import add_residues

def main():
    args = docopt.docopt(__doc__)

    input_path = args['<file_or_folder>']
    resis = args['<residue_string>']
    position = int(args['<pdb_position>'])
    chain = args['--chain']

    init()

    if os.path.isdir(input_path):
        bk_folder = os.path.join(input_path, '..', 'outputs_bk')
        os.mkdir(bk_folder)
        for path in glob.glob(input_path + '/*.pdb*'):
            basename = os.path.basename(path)
            shutil.copyfile(path, os.path.join(bk_folder, basename))
            pose = pose_from_file(path)
            add_residues.add_aas(pose, position, resis, pdbnum=True,
                    chain=chain)
            pose.dump_pdb(path)

    elif os.path.isfile(input_path):
        if not args['--out']:
            shutil.copyfile(input_path, input_path + '.bk')
        pose = pose_from_file(input_path)
        add_residues.add_aas(pose, position, resis, pdbnum=True,
                chain=chain)
        if not args['--out']:
            pose.dump_pdb(input_path)
        else:
            pose.dump_pdb(args['--out'])


