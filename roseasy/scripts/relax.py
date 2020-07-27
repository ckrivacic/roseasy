from pyrosetta import init
from roseasy.workspace import pipeline
from roseasy import big_jobs
import os, sys, subprocess
from roseasy.movers import relax as r

class Workspace(pipeline.RelaxModels):
    pass

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

    pose.dump_pdb(output_prefix + 'input.pdb')
