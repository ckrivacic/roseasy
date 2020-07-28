from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta.core.scoring import CA_rmsd
from pyrosetta.core.scoring import all_atom_rmsd
from roseasy.workspace import pipeline
from roseasy import big_jobs
import os, sys, subprocess
from roseasy.movers import relax as r

def get_workspace(root_dir, step):
    return pipeline.RelaxModels(root_dir, step)

if __name__=='__main__':
    workspace, job_id, task_id, parameters = big_jobs.initiate()
    output_prefix = '{0}/{1}_{2:06d}_'.format(workspace.output_dir,
            job_id, task_id)
    test_run = parameters.get('test_run', False)
    init()
    pdbpath = workspace.input_pdb_path
    pose = pose_from_file(workspace.input_pdb_path)
    relax = r.Relax()
    if test_run:
        relax.rounds = 1
    relax.pose = pose
    relax.apply()

    input_pose = pose_from_file(workspace.input_pdb_path)
    ca_rmsd = CA_rmsd(relax.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(relax.pose, input_pose)

    pose.dump_pdb(output_prefix + 'input.pdb')
    with open(output_prefix + 'input.pdb', 'a') as f:
        f.write('\nCA_RMSD {}'.format(ca_rmsd))
        f.write('\nAllAtom_RMSD {}'.format(all_atom_rmsd))
