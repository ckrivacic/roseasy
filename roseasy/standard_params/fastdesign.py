from pyrosetta import init
from pyrosetta.rosetta.core.pose import setPoseExtraScore
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from roseasy import pipeline
from roseasy import big_jobs
import os, sys, subprocess, gzip
from roseasy.movers import fastdesign
from roseasy.standard_params.filters import FilterContainer

def get_workspace(root_dir, step):
    return pipeline.DesignWorkspace(root_dir, step)

if __name__=='__main__':
    workspace, job_info = big_jobs.initiate()
    test_run = job_info.get('test_run', False)
    init()
    pdbpath = workspace.input_path(job_info)
    pose = pose_from_file(workspace.input_pdb_path)
    fd = fastdesign.FastDesign()

    dalphaball_path = os.path.join(workspace.rosetta_dir, 'source',
            'external', 'DAlphaBall', 'DAlphaBall.gcc')
    fd.add_init_arg('-holes:dalphaball {} -in:file:s {}'.format(dalphaball_path, pdbpath))
    if test_run:
        fd.rounds = 1
    fd.pose = pose
    # Warning: This is an all-atom movemap. Constrain to input coords if
    # you don't want things to move around a lot.
    loop = workspace.largest_loop
    fd.task_factory_from_range(loop.start, loop.end)
    print(fd.movemap)
    print(fd.task_factory)
    fd.apply()

    input_pose = pose_from_file(workspace.input_pdb_path)
    ca_rmsd = CA_rmsd(fd.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(fd.pose, input_pose)
    score_fragments = os.path.exists(workspace.loops_path)

    filters = FilterContainer(workspace, fd.pose,
            task_id=job_info['task_id'], score_fragments=score_fragments,
            test_run=test_run)
    filters.run_filters()

    input_name = os.path.basename(pdbpath).split(".")[0]
    out = workspace.output_prefix(job_info) + input_name + workspace.output_suffix(job_info) + '.pdb.gz'

    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_AllAtom_RMSD', all_atom_rmsd)

    pose.dump_pdb(out)
