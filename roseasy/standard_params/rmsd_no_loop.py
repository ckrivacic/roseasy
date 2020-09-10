'''
Find RMSD of all residues not included in the modelled loop for all files in a folder 
and compare to a native or input structure.

Usage:
    rmsd_no_loop.py <folder> <input>

'''

import docopt
from pyrosetta.rosetta.core.pose import setPoseExtraScore, chain_end_res
from pyrosetta.rosetta.core.scoring import CA_rmsd, native_CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta.rosetta.std import map_unsigned_long_unsigned_long
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

    #mobile = pose_from_file('/Users/benjaminrjagger/UCSF/cas_des/CPP_sims/4un3_1051_penetratin/input.pdb.gz')
    
    # ins_len = chain_end_res(mobile, 1) - chain_end_res(target, 1) 
    # des_res = list(range(1, int(loop.start)-1)) + list(range(int(loop.end)+1, chain_end_res(mobile, 1)))
    # wt_res = list(range(1,int(loop.start) -1)) + list(range(int(loop.end)+1 - ins_len, chain_end_res(target, 1)))

    # res_map = map_unsigned_long_unsigned_long()
    # for i in range(len(des_res)):
    #     res_map[des_res[i]] = wt_res[i]

    # rmsd = CA_rmsd(mobile, target, res_map)
    



    for root, dirs, files in os.walk(workspace.output_dir):
        for name in files:
            if name.endswith('.pdb.gz'):
                pdbpath = os.path.join(root, name)
                mobile = pose_from_file(pdbpath) 
                ins_len = chain_end_res(mobile, 1) - chain_end_res(target, 1) 
                des_res = list(range(1, int(loop.start)-1)) + list(range(int(loop.end)+1, chain_end_res(mobile, 1)))
                wt_res = list(range(1,int(loop.start) -1)) + list(range(int(loop.end)+1 - ins_len, chain_end_res(target, 1)))

                res_map = map_unsigned_long_unsigned_long()
                for i in range(len(des_res)):
                    res_map[des_res[i]] = wt_res[i]

                rmsd = CA_rmsd(mobile, target, res_map)
                metric_name = 'EXTRA_METRIC_CA_RMSD_NO_LOOP [[-]]'

                add_lines_to_gzip(pdbpath, [metric_name + ' ' + str(rmsd)])
                #rmsd = all_atom_rmsd(mobile, target)
                print(rmsd)


if __name__=='__main__':
    main()

