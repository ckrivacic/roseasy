from pyrosetta import init
from pyrosetta.rosetta.core.pose import setPoseExtraScore
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from roseasy import pipeline
from roseasy import big_jobs
import os, sys, subprocess, gzip
from roseasy.movers import relax as r
#from roseasy.standard_params.filters import FilterContainer

def get_workspace(root_dir, step):
    return pipeline.DesignWorkspace(root_dir, step)

if __name__=='__main__':
    workspace, job_info = big_jobs.initiate()
    test_run = job_info.get('test_run', False)
    init('-total_threads 1')
    pdbpath = workspace.input_path(job_info)
    pose = pose_from_file(pdbpath)
    relax = r.Relax()

    dalphaball_path = os.path.join(workspace.rosetta_dir, 'source',
            'external', 'DAlpahBall', 'DAlphaBall.gcc')
    relax.add_init_arg('-holes:dalphaball {} -in:file:s {}'.format(dalphaball_path, pdbpath))
    relax.add_init_arg('-total_threads 1')
    if test_run:
        relax.rounds = 1
    relax.pose = pose
    # Warning: This is an all-atom movemap. Constrain to input coords if
    # you don't want things to move around a lot.
    relax.setup_default_movemap()
    relax.apply()

    # This will compare it to the input to this step
    input_pose = pose_from_file(pdbpath)
    # But you can uncomment this to compare to the input to the
    # workspace
    #input_pose = pose_from_file(workspace.input_pdb_path)
    ca_rmsd = CA_rmsd(relax.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(relax.pose, input_pose)
    score_fragments = os.path.exists(workspace.loops_path)

    filters = workspace.get_filters(relax.pose, 
            task_id=job_info['task_id'], score_fragments=score_fragments,
            test_run=test_run)
    filters.run_filters()

    input_name = os.path.basename(pdbpath).split(".")[0]
    out = workspace.output_prefix(job_info) + input_name + workspace.output_suffix(job_info) + '.pdb.gz'

    setPoseExtraScore(relax.pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
    setPoseExtraScore(relax.pose, 'EXTRA_METRIC_AllAtom_RMSD', all_atom_rmsd)

    pose.dump_pdb(out)
