from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from roseasy import pipeline
from roseasy import big_jobs
import os, sys, subprocess, gzip
from roseasy.movers import relax as r

def get_workspace(root_dir, step):
    return pipeline.DesignWorkspace(root_dir, step)

if __name__=='__main__':
    workspace, job_info = big_jobs.initiate()
    output_prefix = '{0}/{1}_{2:06d}_'.format(workspace.output_dir,
            job_info['job_id'], job_info['task_id'])
    test_run = job_info.get('test_run', False)
    init()
    pdbpath = workspace.input_pdb_path
    pose = pose_from_file(workspace.input_pdb_path)
    relax = r.Relax()
    if test_run:
        relax.rounds = 1
    relax.pose = pose
    relax.setup_default_movemap()
    relax.apply()

    input_pose = pose_from_file(workspace.input_pdb_path)
    ca_rmsd = CA_rmsd(relax.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(relax.pose, input_pose)

    pose.dump_pdb(output_prefix + 'input.pdb.gz')
    with gzip.open(output_prefix + 'input.pdb.gz', 'at') as f:
        f.write('\nEXTRA_METRIC_CA_RMSD {}'.format(ca_rmsd))
        f.write('\nEXTRA_METRIC_AllAtom_RMSD {}'.format(all_atom_rmsd))
